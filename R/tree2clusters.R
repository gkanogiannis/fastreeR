#' Perform Hierarchical Clustering and tree prunning on a phylogenetic tree
#'
#' The phylogenetic tree is pruned with
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} to get clusters.
#'
#' @param tree.str A \code{\link[base]{character} vector} of a
#'   phylogenetic tree in Newick format
#' @param tree.distances code{numeric \link[base]{matrix}} of
#'   distances, that were used to generate the tree.
#'   If NULL, it will be infered from tree branch lengths.
#' @param tree.labels A \code{\link[base]{character} vector} of tree labels.
#' @param cutHeight Define cutHeight for tree cutting.
#'   Default automatically defined.
#' @param minClusterSize Minimum size of clusters. Default 1.
#' @param extra Boolean whether to use extra parameters
#'   for the \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @usage tree2clusters(
#'     tree.str,
#'     tree.distances = NULL,
#'     tree.labels = NULL,
#'     cutHeight = NULL,
#'     minClusterSize = 1,
#'     extra = TRUE
#' )
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
#' my.clust <- dist2clusters(
#'     input = dist2tree(
#'         input = system.file("extdata", "samples.dist.gz",
#'         package = "fastreeR"
#'         )
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}
tree2clusters <- function(
                        tree.str,
                        tree.distances = NULL,
                        tree.labels = NULL,
                        cutHeight = NULL,
                        minClusterSize = 1,
                        extra = TRUE) {
    if (is.null(tree.str)) {
        return(NA)
    }

    hierarchicalcluster <- rJava::J(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )

    # dynamicTreeCut expects 2x the provided cutHeight
    if(!is.null(cutHeight)) {
        cutHeight <- 2.0*cutHeight
    }
    clusters <- dynamicTreeCut(
        tree.str = tree.str,
        distancesReordered = tree.distances,
        labelsReordered = tree.labels,
        minClusterSize = minClusterSize,
        cutHeight = cutHeight,
        extra = extra,
        hierarchicalcluster = hierarchicalcluster
    )
    clusters.vec <- hierarchicalcluster$hclusteringClustersNoJRI(clusters)
    gc()
    rJava::J("java.lang.Runtime")$getRuntime()$gc()
    return(list(tree.str, clusters.vec))
}

dynamicTreeCut <- function(
                        tree.str, distancesReordered = NULL,
                        labelsReordered = NULL, minClusterSize = 1,
                        cutHeight = NULL, extra = FALSE,
                        hierarchicalcluster) {
    if (is.null(tree.str)) {
        return(NA)
    }

    tree <- ape::read.tree(text = tree.str)
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
        treecut_hybrid <- dynamicTreeCut::cutreeDynamic(
            dendro = tree.hclust, cutHeight = cutHeight,
            minClusterSize = minClusterSize, method = "hybrid",
            distM = distancesMatrix, deepSplit = 1, verbose = 2, indent = 0
        )
    } else {
        treecut_hybrid <- dynamicTreeCut::cutreeDynamic(
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
        treecut_hybrid, rownames(distancesMatrix)
    )

    return(clusters)
}
