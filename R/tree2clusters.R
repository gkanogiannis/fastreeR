#' Perform Hierarchical Clustering and tree pruning on a phylogenetic tree
#'
#' The phylogenetic tree is pruned with
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} to get clusters.
#'
#' @param treeStr A \code{\link[base]{character} vector} of a
#'   phylogenetic tree in Newick format
#' @param treeDistances \code{numeric \link[base]{matrix}} of
#'   distances, that were used to generate the tree.
#'   If NULL, it will be inferred from tree branch lengths.
#' @param treeLabels A \code{\link[base]{character} vector}
#'   of tree leaf labels.
#' @param cutHeight Define at which height to cut tree.
#'   Default automatically defined.
#' @param minClusterSize Minimum size of clusters. Default 1.
#' @param extra Boolean whether to use extra parameters
#'   for the \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @return
#' \itemize{
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
#' my.clust <- tree2clusters(
#'     treeStr = dist2tree(
#'         inputDist = system.file("extdata", "samples.vcf.dist.gz",
#'             package = "fastreeR"
#'         )
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

tree2clusters <- function(treeStr, treeDistances = NULL, treeLabels = NULL,
                        cutHeight = NULL, minClusterSize = 1, extra = TRUE) {

    tree2clusters_checkParams(treeStr = treeStr, treeDistances = treeDistances,
                                treeLabels = treeLabels, cutHeight = cutHeight,
                                minClusterSize = minClusterSize, extra = extra)

    hierarchicalcluster <- rJava::.jnew(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )

    # dynamicTreeCut expects 2x the provided cutHeight
    if(!is.null(cutHeight)) {
        cutHeight <- 2.0*cutHeight
    }
    clusters <- dynamicTreeCut(
        treeStr = treeStr,
        distancesReordered = treeDistances,
        labelsReordered = treeLabels,
        minClusterSize = minClusterSize,
        cutHeight = cutHeight,
        extra = extra,
        hierarchicalcluster = hierarchicalcluster
    )
    clusters.vec <- hierarchicalcluster$hclusteringClustersNoJRI(clusters)

    return(clusters.vec)
}

dynamicTreeCut <- function(treeStr, distancesReordered = NULL,
                        labelsReordered = NULL, minClusterSize = 1,
                        cutHeight = NULL, extra = FALSE, hierarchicalcluster) {
    tree <- ape::read.tree(text = treeStr)
    ifelse(ape::is.ultrametric(tree),
        "Is ultrametric=TRUE", "Is ultrametric=FALSE"
    )
    tree.hclust <- stats::as.hclust(tree)

    if (is.null(distancesReordered)) {
        distancesMatrix <- ape::cophenetic.phylo(tree)
    } else {
        distancesMatrix <- distancesReordered
        rawLabels <- labelsReordered
        colLabels <- labelsReordered
        rownames(distancesMatrix) <- rawLabels
        colnames(distancesMatrix) <- colLabels
    }

    if (!extra) {
        treecut.hybrid <- dynamicTreeCut::cutreeDynamic(
            dendro = tree.hclust, cutHeight = cutHeight,
            minClusterSize = minClusterSize, method = "hybrid",
            distM = distancesMatrix, deepSplit = 1, verbose = 2, indent = 0
        )
    } else {
        treecut.hybrid <- dynamicTreeCut::cutreeDynamic(
            dendro = tree.hclust, cutHeight = cutHeight,
            minClusterSize = minClusterSize, method = "hybrid",
            distM = distancesMatrix, deepSplit = 1, maxCoreScatter = NULL,
            minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL,
            minSplitHeight = NULL, minAbsSplitHeight = cutHeight,
            pamStage = TRUE, pamRespectsDendro = TRUE, useMedoids = FALSE,
            maxDistToLabel = NULL, maxPamDist = cutHeight,
            respectSmallClusters = TRUE, verbose = 2, indent = 0
        )
    }

    clusters <- hierarchicalcluster$findClusters(
        treecut.hybrid, rownames(distancesMatrix)
    )

    return(clusters)
}

tree2clusters_checkParams <- function(treeStr, treeDistances, treeLabels,
                                            cutHeight, minClusterSize, extra) {
    if (is.null(treeStr) || !methods::is(treeStr, "character")) {
        stop("treeStr parameter must be a character vector.")
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

    if (!is.null(treeDistances) && !methods::is(treeDistances, "matrix")) {
        stop("treeDistances parameter must be a matrix.")
    }

    if (!is.null(treeLabels) && !methods::is(treeLabels, "character")) {
        stop("treeLabels parameter must be a character vector.")
    }
}
