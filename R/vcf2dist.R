#' Calculate distances between samples of a VCF file
#'
#' This function calculates a cosine type dissimilarity measurement between the
#' \code{n} samples of a VCF file.
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
#' Output file, if provided, will contain \code{n+1} lines.
#' The first line contains the number \code{n} of samples
#' and number \code{m} of variants, separated by space.
#' Each of the subsequent \code{n} lines contains \code{n+1} values,
#' separated by space.
#' The first value of each line is a sample name
#' and the rest \code{n} values
#' are the calculated distances of this sample to all the samples.
#' Example output file of the distances of 3 samples
#' calculated from 1000 variants:
#' \tabular{llll}{
#'     3 1000 \tab \cr
#'     Sample1 \tab 0.0 \tab 0.5 \tab 0.2\cr
#'     Sample2 \tab 0.5 \tab 0.0 \tab 0.9\cr
#'     Sample3 \tab 0.2 \tab 0.9 \tab 0.0\cr
#' }
#'
#' @param inputfile Input vcf file location (uncompressed or gzip compressed).
#' @param outputfile Output distances file location.
#' @param threads Number of java threads to use.
#' @param ignoremissing Ignore variants with missing data
#'     (\code{./.} or \code{.|.})
#' @param onlyhets Only calculate on variants with heterozygous calls.
#' @param ignorehets Only calculate on variants with homozygous calls.
#' @param compress Compress output (adds .gz extension).
#'
#' @usage
#' vcf2dist(
#'     inputfile,
#'     outputfile = NULL,
#'     threads = 1,
#'     ignoremissing = FALSE,
#'     onlyhets = FALSE,
#'     ignorehets = FALSE,
#'     compress = TRUE
#' )
#'
#' @return A \code{\link[stats]{dist}} distances object of the calculation.
#' @export
#'
#' @examples
#' my.dist <- vcf2dist(
#'     inputfile = system.file("extdata", "samples.vcf.gz",
#'         package = "fastreeR"
#'     ),
#'     threads = 1
#' )
#' plot(stats::hclust(my.dist))
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

vcf2dist <- function(inputfile, outputfile = NULL, threads = 1,
                    ignoremissing = FALSE, onlyhets = FALSE,
                    ignorehets = FALSE, compress = TRUE) {
    if (is.null(inputfile) || !file.exists(inputfile)) {return(NA)}

    if (R.utils::isGzipped(inputfile)) {
        temp.in <- tempfile(fileext = ".vcf")
        on.exit(unlink(temp.in))
        R.utils::gunzip(
            filename = inputfile, destname = temp.in, remove = FALSE)
        inputfile <- temp.in
    }

    bioinfojavautils <- rJava::.jnew(
        "ciat/agrobio/javautils/JavaUtils", class.loader = .rJava.class.loader)

    cmd <- paste("VCF2DIST", "--numberOfThreads", threads,
            ifelse(ignoremissing, "--ignoremissing", ""),
            ifelse(onlyhets, "--onlyHets", ""),
            ifelse(ignorehets, "--ignoreHets", ""), inputfile, sep = " ")

    temp.out <- tempfile(fileext = ".txt")
    on.exit(unlink(temp.out))
    jSys <- rJava::J("java/lang/System")
    jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$main(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")
    ret.df <- utils::read.table(text = ret.str[-1])
    ret.names <- ret.df[, 1]
    ret.df <- ret.df[, -1]
    rownames(ret.df) <- ret.names
    colnames(ret.df) <- ret.names

    if (!is.null(outputfile)) {
        if (compress) {
            temp.dist <- tempfile(fileext = ".dist")
            on.exit(unlink(temp.dist))
            data.table::fwrite(as.list(ret.str), file = temp.dist, sep = "\n")
            R.utils::gzip(filename = temp.dist,
                destname = paste0(outputfile, ".gz"), overwrite = TRUE
            )
        } else {
            data.table::fwrite(as.list(ret.str), file = outputfile, sep = "\n")
        }
    }
    return(stats::as.dist(as.matrix(ret.df)))
}
