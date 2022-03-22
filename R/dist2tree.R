#' Generate phylogenetic tree from samples of a distance matrix
#'
#' Performs Hierarchical Clustering on a distance matrix
#' (i.e. calculated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}})
#' and generates a phylogenetic tree with
#' agglomerative Neighbor Joining method (complete linkage).
#'
#' @param input.dist Input distances file location
#' (generated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}}).
#' File can be gzip compressed.
#' Or a \code{\link[stats]{dist}} distances object.
#'
#' @usage
#' dist2tree(
#'     input.dist
#' )
#'
#' @return A \code{\link[base]{character} vector} of the generated
#' phylogenetic tree in Newick format.
#' @export
#'
#' @examples
#' my.tree <- dist2tree(
#'     input.dist =
#'         system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2tree <- function(input.dist) {
    if (is.null(input.dist) ||
        (!methods::is(input.dist, "dist") &&
            !methods::is(input.dist, "character")) ||
        (methods::is(input.dist, "character") && !file.exists(input.dist))) {
        invisible(NULL)
    }

    inputfile <- input.dist

    if(methods::is(input.dist, "character") && R.utils::isGzipped(input.dist)) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        R.utils::gunzip(filename = input.dist, destname = temp.in,
                        remove = FALSE)
        inputfile <- temp.in
    } else if (methods::is(input.dist, "dist")) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        write(paste0(attr(input.dist, "Size"), " 0"), file = temp.in)
        utils::write.table(
            as.matrix(input.dist),
            file = temp.in,
            sep = " ",
            row.names = TRUE,
            col.names = FALSE,
            append = TRUE
        )
        inputfile <- temp.in
    }

    hierarchicalcluster <- rJava::J(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )
    generaltools <- rJava::J(
        class="ciat/agrobio/core/GeneralTools",
        class.loader = .rJava.class.loader
    )

    data <- generaltools$readDistancesSamples(inputfile)
    samples.names <- rJava::.jevalArray(data[[2]])
    samples.distances <- rJava::.jevalArray(data[[1]], simplify = TRUE)

    tree.str <- hierarchicalcluster$hclusteringTree(data[[2]], data[[1]])

    return(tree.str)
}
