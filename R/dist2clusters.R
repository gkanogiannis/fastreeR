#' Perform Hierarchical Clustering and tree prunning on a distance matrix
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
#' @param input.dist Input distances file location
#' (generated with \code{\link[fastreeR]{vcf2dist}}
#' or \code{\link[fastreeR]{fasta2dist}}).
#' File can be gzip compressed.
#' Or a \code{\link[stats]{dist}} distances object.
#' @param cutHeight Define cutHeight for tree cutting.
#' Default automatically defined.
#' @param minClusterSize Minimum size of clusters. Default 1.
#' @param extra Boolean whether to use extra parameters
#' for the \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @usage
#' dist2clusters(
#'     input.dist,
#'     cutHeight = NULL,
#'     minClusterSize = 1,
#'     extra = TRUE
#' )
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
#'     input.dist =
#'         system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2clusters <- function(input.dist, cutHeight = NULL, minClusterSize = 1,
                            extra = TRUE) {
    if (is.null(input.dist) ||
        (!methods::is(input.dist, "dist") &&
            !methods::is(input.dist, "character")) ||
        (methods::is(input.dist, "character") &&
            (!file.exists(input.dist) || length(input.dist==0)))) {
        return(NULL)
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
        utils::write.table(as.matrix(input.dist), file = temp.in, sep = " ",
                            row.names = TRUE, col.names = FALSE, append = TRUE)
        inputfile <- temp.in
    }

    hierarchicalcluster <- rJava::J(
        class="ciat/agrobio/hcluster/HierarchicalCluster",
        class.loader = .rJava.class.loader
    )
    generaltools <- rJava::J(
        class="ciat/agrobio/core/GeneralTools",
        class.loader = .rJava.class.loader)

    # data[[1]] distances, data[[2]] labels
    data <- generaltools$readDistancesSamples(inputfile)
    tree.str <- hierarchicalcluster$hclusteringTree(data[[2]], data[[1]])

    labelsReordered <- generaltools$reorderLabels(data[[2]], tree.str)
    distancesReordered <- rJava::.jevalArray(
        generaltools$reorderDistances(data[[1]], data[[2]], labelsReordered),
        simplify = TRUE
    )

    return(list(
        tree.str,
        fastreeR::tree2clusters(
            tree.str = tree.str, tree.distances = distancesReordered,
            tree.labels = labelsReordered, cutHeight = cutHeight,
            minClusterSize = minClusterSize, extra = extra
        )
    ))
}
