#' Calculate various per sample statistics from a VCF file
#'
#' Only biallelic SNPs are considered.
#' For each sample, the following statistics are calculated :
#' \itemize{
#'     \item INDIV : Sample name
#'     \item N_SITES : Total number of SNPs
#'     \item N_HET : Number of SNPs with
#'         heterozygous call (\code{0/1} or \code{0|1}
#'                         or \code{1/0} or \code{1|0})
#'     \item N_ALT : Number of SNPs with
#'         alternate homozygous call (\code{1/1} or \code{1|1})
#'     \item N_REF : Number of SNPs with
#'         reference homozygous call (\code{0/0} or \code{0|0})
#'     \item N_MISS : Number of SNPs with
#'         missing call (\code{./.} or \code{.|.})
#'     \item P_HET : Percentage of heterozygous calls
#'     \item P_ALT : Percentage of alternate homozygous calls
#'     \item P_REF : Percentage of reference homozygous calls
#'     \item P_MISS : Percentage of missing calls (missing rate)
#' }
#'
#' @param inputFile Input vcf file location (uncompressed or gzip compressed).
#' @param outputFile Output samples statistics file location.
#'
#' @return A \code{\link[base]{data.frame}} of sample statistics.
#' @export
#'
#' @examples
#' my.istats <- vcf2istats(
#'     inputFile =
#'         system.file("extdata", "samples.vcf.gz", package = "fastreeR")
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

vcf2istats <- function(inputFile, outputFile = NULL) {
    vcf2istats_checkParams(inputFile = inputFile, outputFile = outputFile)

    bioinfojavautils <- rJava::.jnew(
        class="ciat/agrobio/javautils/JavaUtils",
        class.loader = .rJava.class.loader
    )
    cmd <- paste(
        "VCF2ISTATS",
        inputFile,
        sep = " "
    )

    temp.out <- tempfile(fileext = ".txt")
    on.exit(unlink(temp.out))
    jSys <- rJava::J("java/lang/System")
    jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")
    ret.df <- utils::read.table(text = ret.str, header = TRUE)

    if (!is.null(outputFile)) {
        data.table::fwrite(as.list(ret.str), file = outputFile, sep = "\n")
    }

    return(ret.df)
}

vcf2istats_checkParams <- function(inputFile, outputFile) {
    if (!methods::is(inputFile, "character")){
        stop("inputFile must be a file location.")
    }

    if (is.null(inputFile) || !file.exists(inputFile)) {
        stop("inputFile=",inputFile," does not exist.")
    }

    if ((!is.null(outputFile) && !methods::is(outputFile, "character")) ||
        (methods::is(outputFile, "character") && nchar(outputFile)==0)) {
        stop("outputFile must be a file location.")
    }
}
