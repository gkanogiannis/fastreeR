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
#' @param verbose Logical. If TRUE, enables verbose output from the Java backend.
#'
#' @return A \code{\link[base]{character} vector} of the generated
#' phylogenetic tree in Newick format.
#' @export
#'
#' @examples
#' my.tree <- dist2tree(
#'     inputDist =
#'     system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR"),
#'     verbose = TRUE
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2tree <- function(inputDist, verbose = FALSE) {
    dist2tree_checkParams(inputDist = inputDist, verbose = verbose)

    inputfile <- inputDist

    if (methods::is(inputDist, "character") && R.utils::isGzipped(inputDist)) {
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
        class="com/gkano/bioinfo/tree/HierarchicalCluster",
        verbose,
        class.loader = .rJava.class.loader
    )
    generaltools <- rJava::J(
        class="com/gkano/bioinfo/var/GeneralTools",
        class.loader = .rJava.class.loader
    )$getInstance()

    # data[[1]] distances, data[[2]] labels
    data <- generaltools$readDistancesSamples(inputfile)
    
    # Call the Java method which returns an Object[] and convert it to R
    ret_java <- hierarchicalcluster$hclusteringTree(data[[2]], data[[1]])
    # Convert Java array to an R vector; simplify=TRUE will convert Java strings to R character
    ret_r <- rJava::.jevalArray(ret_java, simplify = TRUE)
    if (length(ret_r) >= 1 && is.character(ret_r[[1]])) {
        # ret_r is a character vector when Java returned strings
        treeStr <- as.character(ret_r[[1]])
    } else if (length(ret_r) >= 1) {
        # Fallback: try calling Java toString() on the first element
        first_java_elem <- ret_java[[1]]
        # safe try
        treeStr <- tryCatch(
            rJava::.jcall(first_java_elem, "S", "toString"),
            error = function(e) character(0)
        )
    } else {
        treeStr <- character(0)
    }

    return(treeStr)
}

dist2tree_checkParams <- function(inputDist, verbose) {
    if (is.null(inputDist) ||
        (!methods::is(inputDist, "dist") &&
         !methods::is(inputDist, "character")) ||
        (methods::is(inputDist, "character") &&
         (!file.exists(inputDist) || nchar(inputDist)==0))) {
        stop("inputDist parameter must be a valid file location ",
             "or a dist object.")
    }

    if (!is.logical(verbose)){
        stop("verbose",
             "must be logical.")
    }
}
