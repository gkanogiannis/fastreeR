#' Generate phylogenetic tree from samples of a VCF file
#'
#' This function calculates a distance matrix between the samples of a VCF file
#' as in \code{\link[fastreeR]{vcf2dist}}
#' and performs Hierarchical Clustering on this distance matrix
#' as in \code{\link[fastreeR]{dist2tree}}.
#' A phylogenetic tree is calculated with
#' agglomerative Neighbor Joining method (complete linkage).
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
#'
#' @usage
#' vcf2tree(
#'     inputfile,
#'     threads = 1,
#'     ignoremissing = FALSE,
#'     onlyhets = FALSE,
#'     ignorehets = FALSE
#' )
#'
#' @return A \code{\link[base]{character}} vector of the generated
#' phylogenetic tree in Newick format.
#' @export
#'
#' @examples
#' my.tree <- vcf2tree(
#'     inputfile = system.file("extdata", "samples.vcf.gz",
#'         package = "fastreeR"
#'     ),
#'     threads = 1
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

vcf2tree <- function(
                    inputfile,
                    threads = 1,
                    ignoremissing = FALSE,
                    onlyhets = FALSE,
                    ignorehets = FALSE) {
    if (is.null(inputfile) || !file.exists(inputfile)) {
        return(NA)
    }

    if (R.utils::isGzipped(inputfile)) {
        temp.in <- tempfile(fileext = ".vcf")
        on.exit(unlink(temp.in))
        R.utils::gunzip(filename = inputfile,
                        destname = temp.in,
                        remove = FALSE)
        inputfile <- temp.in
    }

    bioinfojavautils <- rJava::.jnew(
        "ciat/agrobio/javautils/JavaUtils",
        class.loader = .rJava.class.loader
    )

    cmd <- paste(
        "VCF2TREE",
        "--numberOfThreads", threads,
        ifelse(ignoremissing, "--ignoremissing", ""),
        ifelse(onlyhets, "--onlyHets", ""),
        ifelse(ignorehets, "--ignoreHets", ""),
        inputfile,
        sep = " "
    )

    temp.out <- tempfile(fileext = ".txt")
    on.exit(unlink(temp.out))
    jSys <- rJava::J("java/lang/System")
    jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$main(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")

    return(ret.str)
}
