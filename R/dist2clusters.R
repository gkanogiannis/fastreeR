#' Perform Hierarchical Clustering and tree pruning on a distance matrix
#'
#' Performs Hierarchical Clustering on a distance matrix
#' (i.e. calculated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}})
#' and generates a phylogenetic tree with
#' agglomerative Neighbor Joining method (complete linkage)
#' (as in \code{\link[fastreeR]{dist2tree}}).
#' The phylogenetic tree is then pruned with
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} to get clusters
#' (as in \code{\link[fastreeR]{tree2clusters}}).
#'
#' @param inputDist Input distances file location
#' (generated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}}).
#' File can be gzip compressed.
#' Or a \code{\link[stats]{dist}} distances object.
#' @param cutHeight Define at which height to cut tree.
#' Default automatically defined.
#' @param minClusterSize Minimum size of clusters. Default 1.
#' @param extra Boolean whether to use extra parameters
#' for the \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @return A list of :
#' \itemize{
#'     \item \code{\link[base]{character} vector} of the generated
#'     phylogenetic tree in Newick format
#'     \item \code{\link[base]{character} vector} of the clusters.
#'     Each row contains data for a cluster, separated by space.
#'     The id of the cluster,
#'     the size of the cluster (number of elements)
#'     and the names of its elements,
#'     Cluster id 0 contains all the objects not assigned
#'     to a cluster (singletons).
#'     Example clusters output :
#'     \tabular{lllll}{
#'         0 \tab 3 \tab Sample1 \tab Sample2 \tab Sample3 \cr
#'         1 \tab 3 \tab Sample4 \tab Sample5 \tab Sample6 \cr
#'         2 \tab 2 \tab Sample7 \tab Sample8 \tab \cr
#'         3 \tab 2 \tab Sample9 \tab Sample0\tab \cr
#'     }
#' }
#' @export
#'
#' @examples
#' my.clust <- dist2clusters(
#'     inputDist =
#'         system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2clusters <- function(inputDist, cutHeight = NULL,
                                            minClusterSize = 1, extra = TRUE){
    dist2clusters_checkParams(inputDist = inputDist, cutHeight = cutHeight,
                                minClusterSize = minClusterSize, extra = extra)

    inputfile <- inputDist

    if(methods::is(inputDist, "character") && R.utils::isGzipped(inputDist)) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        R.utils::gunzip(filename = inputDist, destname = temp.in,remove = FALSE)
        inputfile <- temp.in
    } else if (methods::is(inputDist, "dist")) {
        temp.in <- tempfile(fileext = ".dist")
        on.exit(unlink(temp.in))
        write(paste0(attr(inputDist, "Size"), " 0"), file = temp.in)
        utils::write.table(as.matrix(inputDist), file = temp.in, sep = " ",
            quote = FALSE, row.names = TRUE, col.names = FALSE, append = TRUE)
        inputfile <- temp.in
    }

    hierarchicalcluster <- rJava::.jnew(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )
    generaltools <- rJava::J(class="ciat/agrobio/core/GeneralTools",
        class.loader = .rJava.class.loader)$getInstance()

    # data[[1]] distances, data[[2]] labels
    data <- generaltools$readDistancesSamples(inputfile)
    treeStr <- hierarchicalcluster$hclusteringTree(data[[2]], data[[1]])

    labelsReordered <- generaltools$reorderLabels(data[[2]], treeStr)
    distancesReordered <- rJava::.jevalArray(
        generaltools$reorderDistances(data[[1]], data[[2]], labelsReordered),
        simplify = TRUE
    )
    return(list(treeStr, fastreeR::tree2clusters(
                        treeStr = treeStr, treeDistances = distancesReordered,
                        treeLabels = labelsReordered, cutHeight = cutHeight,
                        minClusterSize = minClusterSize, extra = extra)))
}

dist2clusters_checkParams <- function(inputDist, cutHeight,
                                                        minClusterSize,extra){
    if (is.null(inputDist) || (!methods::is(inputDist, "dist") &&
            !methods::is(inputDist, "character")) ||
        (methods::is(inputDist, "character") &&
            (!file.exists(inputDist) || nchar(inputDist)==0))) {
        stop("inputDist parameter must be a valid file location ",
                                                            "or a dist object.")
    }

    if (!is.null(cutHeight) &&
        (!is.numeric(cutHeight) || (is.numeric(cutHeight) && cutHeight<0))) {
        stop("cutHeight parameter must be positive numeric.")
    }

    if ((!is.numeric(minClusterSize) ||
        (is.numeric(minClusterSize) && minClusterSize<1))) {
        stop("threads and minClusterSize parameters must be positive integer.")
    }

    if (!is.logical(extra)){
        stop("extra parameter must be logical.")
    }
}
