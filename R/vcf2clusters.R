#' Perform Hierarchical Clustering and tree pruning on samples of VCF file
#'
#' Performs Hierarchical Clustering on a distance matrix
#' calculated as in \code{\link[fastreeR]{vcf2dist}}
#' and generates a phylogenetic tree with
#' agglomerative Neighbor Joining method (complete linkage)
#' (as in \code{\link[fastreeR]{dist2tree}}).
#' The phylogenetic tree is then pruned with
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} to get clusters
#' (as in \code{\link[fastreeR]{tree2clusters}}).
#'
#' Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL
#' variants are considered, phased or not. Some VCF encoding examples are:
#'
#'     \itemize{
#'         \item heterozygous variants : \code{1/0} or \code{0/1} or \code{0/2}
#'         or \code{1|0} or \code{0|1} or \code{0|2}
#'         \item homozygous to the reference allele variants : \code{0/0}
#'         or \code{0|0}
#'         \item homozygous to the first alternate allele variants : \code{1/1}
#'         or \code{1|1}
#'     }
#'
#' If there are \code{n} samples and \code{m} variants, an \code{nxn}
#' zero-diagonal symmetric distance matrix is calculated.
#' The calculated cosine type distance (1-cosine_similarity)/2 is in the range
#' [0,1] where value 0 means completely identical samples (cosine is 1),
#' value 0.5 means perpendicular samples (cosine is 0)
#' and value 1 means completely opposite samples (cosine is -1).
#'
#' The calculation is performed by a Java backend implementation,
#' that supports multi-core CPU utilization
#' and can be demanding in terms of memory resources.
#' By default a JVM is launched with a maximum memory allocation of 512 MB.
#' When this amount is not sufficient,
#' the user needs to reserve additional memory resources,
#' before loading the package,
#' by updating the value of the \code{java.parameters} option.
#' For example in order to allocate 4GB of RAM,
#' the user needs to issue \code{options(java.parameters="-Xmx4g")}
#' before \code{library(fastreeR)}.
#'
#' @param inputfile Input vcf file location (uncompressed or gzip compressed).
#' @param threads Number of java threads to use.
#' @param ignoremissing Ignore variants with missing data
#'     (\code{./.} or \code{.|.})
#' @param onlyhets Only calculate on variants with heterozygous calls.
#' @param ignorehets Only calculate on variants with homozygous calls.
#' @param cutHeight Define cutHeight for tree cutting.
#'     Default automatically defined.
#' @param minClusterSize Minimum size of clusters. Default 1.
#' @param extra Boolean whether to use extra parameters
#'     for the \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @usage
#' vcf2clusters(
#'     inputfile,
#'     threads = 2,
#'     ignoremissing = FALSE,
#'     onlyhets = FALSE,
#'     ignorehets = FALSE,
#'     cutHeight = NULL,
#'     minClusterSize = 1,
#'     extra = TRUE
#' )
#'
#' @return A list of :
#' \itemize{
#'     \item \code{\link[stats]{dist}} distances object.
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
#' my.clust <- vcf2clusters(
#'     inputfile = system.file("extdata", "samples.vcf.gz",
#'         package = "fastreeR"
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}
vcf2clusters <- function(inputfile,
                        threads = 2,
                        ignoremissing = FALSE,
                        onlyhets = FALSE,
                        ignorehets = FALSE,
                        cutHeight = NULL,
                        minClusterSize = 1,
                        extra = TRUE) {
    if (is.null(inputfile) || !file.exists(inputfile)) {
        return(NULL)
    }
    if (R.utils::isGzipped(inputfile)) {
        temp.in <- tempfile(fileext = ".vcf")
        on.exit(unlink(temp.in))
        R.utils::gunzip(
            filename = inputfile,
            destname = temp.in, remove = FALSE
        )
        inputfile <- temp.in
    }

    my.dist <- fastreeR::vcf2dist(
        inputfile = inputfile,
        threads = threads,
        ignoremissing = ignoremissing,
        onlyhets = onlyhets,
        ignorehets = ignorehets
    )
    my.clusters <- fastreeR::dist2clusters(
        input = my.dist,
        cutHeight = cutHeight,
        minClusterSize = minClusterSize,
        extra = extra
    )

    return(list(
        my.dist,
        my.clusters[[1]],
        my.clusters[[2]]
    ))
}
