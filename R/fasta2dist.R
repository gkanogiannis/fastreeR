#' Calculate distances between sequences of a FASTA file
#'
#' This function calculates a d2_S type dissimilarity measurement between the
#' \code{n} sequences (which can represent samples) of a FASTA file.
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input fasta files locations (uncompressed or gzip compressed).
#' @param outputfile Output distances file location.
#' @param threads Number of java threads to use.
#' @param kmer Kmer length to use for analyzing fasta sequences.
#' @param normalize Normalize on sequences length.
#' @param compress Compress output (adds .gz extension).
#'
#' @usage
#' fasta2dist(
#'     ...,
#'     outputfile = NULL,
#'     threads = 2,
#'     kmer = 6,
#'     normalize = FALSE,
#'     compress = TRUE
#' )
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
#' if(!is.null(my.dist))
#'     plot(stats::hclust(my.dist))
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

fasta2dist <- function(...,
                        outputfile = NULL,
                        threads = 2,
                        kmer = 6,
                        normalize = FALSE,
                        compress = TRUE) {
    ins <- unlist(list(...))
    if (length(ins)==0 || list(NULL) %in% ins) {return(NULL)}
    inputfile <- tempfile(fileext = ".fasta")
    on.exit(unlink(inputfile))
    for (i in ins){
        if (!file.exists(i)) {return(NULL)}
        if (R.utils::isGzipped(i)) {
            temp.in <- tempfile(fileext = ".fasta")
            on.exit(unlink(temp.in))
            R.utils::gunzip(filename = i, destname = temp.in, remove = FALSE)
            i <- temp.in
        }
        write(readLines(i), file = inputfile, append = TRUE)
    }

    bioinfojavautils <- rJava::J(
        class="ciat/agrobio/javautils/JavaUtils",
        class.loader = .rJava.class.loader)
    cmd <- paste("FASTA2DIST", "--numberOfThreads", threads,
                ifelse(normalize, "--normalize", ""),
                "--kmerSize", kmer,
                "--inputFile", inputfile, sep = " ")

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
