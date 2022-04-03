#' Generate phylogenetic tree from samples of a distance matrix
#'
#' Performs Hierarchical Clustering on a distance matrix
#' (i.e. calculated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}})
#' and generates a phylogenetic tree with
#' agglomerative Neighbor Joining method (complete linkage).
#'
#' @param inputDist Input distances file location
#' (generated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}}).
#' File can be gzip compressed.
#' Or a \code{\link[stats]{dist}} distances object.
#'
#' @return A \code{\link[base]{character} vector} of the generated
#' phylogenetic tree in Newick format.
#' @export
#'
#' @examples
#' my.tree <- dist2tree(
#'     inputDist =
#'         system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2tree <- function(inputDist) {
    dist2tree_checkParams(inputDist = inputDist)

    inputfile <- inputDist

    if(methods::is(inputDist, "character") && R.utils::isGzipped(inputDist)) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        R.utils::gunzip(filename = inputDist, destname = temp.in,
                        remove = FALSE)
        inputfile <- temp.in
    } else if (methods::is(inputDist, "dist")) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        write(paste0(attr(inputDist, "Size"), " 0"), file = temp.in)
        utils::write.table(
            as.matrix(inputDist),
            file = temp.in,
            sep = " ",
            row.names = TRUE,
            col.names = FALSE,
            append = TRUE
        )
        inputfile <- temp.in
    }

    hierarchicalcluster <- rJava::.jnew(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )
    generaltools <- rJava::J(
        class="ciat/agrobio/core/GeneralTools",
        class.loader = .rJava.class.loader
    )$getInstance()

    data <- generaltools$readDistancesSamples(inputfile)
    samples.names <- rJava::.jevalArray(data[[2]])
    samples.distances <- rJava::.jevalArray(data[[1]], simplify = TRUE)

    treeStr <- hierarchicalcluster$hclusteringTree(data[[2]], data[[1]])

    return(treeStr)
}

dist2tree_checkParams <- function(inputDist) {
    if (is.null(inputDist) ||
        (!methods::is(inputDist, "dist") &&
            !methods::is(inputDist, "character")) ||
        (methods::is(inputDist, "character") &&
            (!file.exists(inputDist) || nchar(inputDist)==0))) {
        stop("inputDist parameter must be a valid file location ",
                                                            "or a dist object.")
    }
}
