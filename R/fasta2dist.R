#' Calculate distances between sequences of a FASTA file
#'
#' This function calculates a d2_S type dissimilarity measurement between the
#' \code{n} sequences (which can represent samples) of a FASTA file.
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input fasta files locations (uncompressed or gzip compressed).
#' @param outputFile Output distances file location.
#' @param threads Number of java threads to use.
#' @param kmer Kmer length to use for analyzing fasta sequences.
#' @param normalize Normalize on sequences length.
#' @param compress Compress output (adds .gz extension).
#'
#' @return A \code{\link[stats]{dist}} distances object of the calculation.
#' @export
#'
#' @examples
#' my.dist <- fasta2dist(
#'     inputfile = system.file("extdata", "samples.fasta.gz",
#'         package = "fastreeR"
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}
#'

fasta2dist <- function(..., outputFile = NULL, threads = 2, kmer = 6,
                        normalize = FALSE, compress = TRUE) {
    ins <- unlist(list(...))

    fasta2dist_checkParams(ins = ins, outputFile = outputFile,
                                threads = threads, kmer = kmer,
                                    normalize = normalize, compress = compress)
    inputfile <- tempfile(fileext = ".fasta")
    on.exit(unlink(inputfile))

    for (i in ins){
        if (!methods::is(i, "character") || is.null(i) || !file.exists(i) ) {
            stop("Input fasta file ",i, " is not valid.")}
        if (R.utils::isGzipped(i)) {
            temp.in <- tempfile(fileext = ".fasta"); on.exit(unlink(temp.in))
            R.utils::gunzip(filename = i, destname = temp.in, remove = FALSE)
            i <- temp.in
        }
        write(readLines(i), file = inputfile, append = TRUE)
    }

    bioinfojavautils <- rJava::.jnew(
        class="ciat/agrobio/javautils/JavaUtils",
        class.loader = .rJava.class.loader)
    cmd <- paste("FASTA2DIST", "--numberOfThreads", threads,
                ifelse(normalize, "--normalize", ""),
                "--kmerSize", kmer, "--inputFile", inputfile, sep = " ")
    temp.out <- tempfile(fileext = ".txt"); on.exit(unlink(temp.out))
    jSys <- rJava::J("java/lang/System"); jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")
    ret.df <- utils::read.table(text = ret.str[-1])
    ret.names <- ret.df[, 1]; ret.df <- ret.df[, -1]
    rownames(ret.df) <- ret.names; colnames(ret.df) <- ret.names
    if (!is.null(outputFile)) {
        if (compress) {
            temp.dist <- tempfile(fileext = ".dist"); on.exit(unlink(temp.dist))
            data.table::fwrite(as.list(ret.str), file = temp.dist, sep = "\n")
            R.utils::gzip(filename = temp.dist,
                        destname = paste0(outputFile, ".gz"), overwrite = TRUE)
        } else {
            data.table::fwrite(as.list(ret.str), file = outputFile, sep = "\n")
        }
    }

    return(stats::as.dist(as.matrix(ret.df), diag = TRUE, upper = TRUE))
}

fasta2dist_checkParams <- function(ins, outputFile, threads, kmer,
                                                        normalize, compress) {
    if (length(ins)==0 || list(NULL) %in% ins) {
        stop("No input fasta files were provided.")
    }

    if ((!is.null(outputFile) && !methods::is(outputFile, "character")) ||
        (methods::is(outputFile, "character") && nchar(outputFile)==0)) {
        stop("outputFile must be a file location.")
    }

    if (!is.numeric(threads) || (is.numeric(threads) && threads<1)) {
        stop("threads parameter must be positive integer.")
    }

    if (!is.numeric(kmer) || (is.numeric(kmer) && kmer<2)) {
        stop("kmer parameter must be positive integer >1 .")
    }

    if(!is.logical(normalize) || !is.logical(compress)) {
        stop("normalize and compress parameters must be logical.")
    }
}
