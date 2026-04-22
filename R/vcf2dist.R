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
#' \strong{Windowed mode.} Setting \code{windowBp} or \code{windowVariants}
#' (mutually exclusive) instructs the Java backend to emit one distance matrix
#' per genomic window. Windows never straddle chromosomes. The return value
#' changes accordingly:
#' \itemize{
#'   \item \code{longFormat = FALSE} (default): a named \code{list} of
#'     \code{\link[stats]{dist}} objects, one per window. List names follow
#'     the format \code{"chrom:start-end"}.
#'   \item \code{longFormat = TRUE}: a single \code{data.frame} with columns
#'     \code{chrom, start, end, sample_i, sample_j, dist}.
#' }
#'
#' @param inputFile Input vcf file location (uncompressed or gzip compressed).
#' @param outputFile Output distances file location.
#' @param threads Number of java threads to use.
#' @param compress Compress output (adds .gz extension).
#' @param verbose Logical. If TRUE, enables verbose output from the Java backend.
#' @param windowBp Optional positive integer. Emit one matrix per window of
#'   N base pairs. Mutually exclusive with \code{windowVariants}.
#' @param windowVariants Optional positive integer. Emit one matrix per N
#'   consecutive variants. Mutually exclusive with \code{windowBp}.
#' @param windowStep Optional positive integer. Window step (defaults to
#'   window size, i.e. tiled). Sliding windows are not yet implemented and
#'   will be rejected by the Java backend.
#' @param windowMinVariants Minimum number of variants required to emit a
#'   window (default 1; smaller windows are skipped silently).
#' @param longFormat Logical. In windowed mode, return a long-form
#'   \code{data.frame} instead of a list of \code{dist} objects.
#'
#' @return In non-windowed mode, a \code{\link[stats]{dist}} distances object.
#'   In windowed mode, either a named \code{list} of \code{dist} objects
#'   (default) or a long-form \code{data.frame} (when \code{longFormat = TRUE}).
#' @export
#'
#' @examples
#' my.dist <- vcf2dist(
#'     inputFile = system.file("extdata", "samples.vcf.gz",
#'         package = "fastreeR"
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

vcf2dist <- function(inputFile, outputFile=NULL,
                    threads=2, compress = FALSE,
                    verbose = FALSE,
                    windowBp = NULL, windowVariants = NULL,
                    windowStep = NULL, windowMinVariants = 1L,
                    longFormat = FALSE) {

    vcf2dist_checkParams(inputFile = inputFile, outputFile = outputFile,
        threads = threads,
        compress = compress,
        verbose = verbose,
        windowBp = windowBp, windowVariants = windowVariants,
        windowStep = windowStep, windowMinVariants = windowMinVariants,
        longFormat = longFormat)

    if (R.utils::isGzipped(inputFile)) {
        temp.in <- tempfile(fileext = ".vcf"); on.exit(unlink(temp.in))
        R.utils::gunzip(filename = inputFile, destname = temp.in, remove=FALSE)
        inputFile <- temp.in
    }

    windowed <- !is.null(windowBp) || !is.null(windowVariants)

    bioinfojavautils <- rJava::.jnew(class="com/gkano/bioinfo/javautils/JavaUtils",
                                    class.loader = .rJava.class.loader)
    cmd_parts <- c("VCF2DIST",
                   "--numberOfThreads", threads,
                   if (verbose) "--verbose" else "",
                   "--input", inputFile)
    if (!is.null(windowBp)) {
        cmd_parts <- c(cmd_parts, "--window-bp", windowBp)
    }
    if (!is.null(windowVariants)) {
        cmd_parts <- c(cmd_parts, "--window-variants", windowVariants)
    }
    if (!is.null(windowStep)) {
        cmd_parts <- c(cmd_parts, "--step", windowStep)
    }
    if (windowed) {
        cmd_parts <- c(cmd_parts, "--min-variants", windowMinVariants)
        if (longFormat) cmd_parts <- c(cmd_parts, "--long")
    }
    cmd <- paste(cmd_parts, collapse = " ")

    temp.out <- tempfile(fileext = ".txt"); on.exit(unlink(temp.out), add = TRUE)
    jSys <- rJava::J("java/lang/System"); jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    if (windowed) {
        raw_lines <- readLines(temp.out)
        if (longFormat) {
            ret_obj <- vcf2dist_parseWindowedLong(raw_lines)
        } else {
            ret_obj <- vcf2dist_parseWindowedConcat(raw_lines)
        }
        if (!is.null(outputFile)) {
            if (compress) {
                temp.dist <- tempfile(fileext = ".dist")
                on.exit(unlink(temp.dist), add = TRUE)
                data.table::fwrite(as.list(raw_lines), file = temp.dist,
                                   sep = "\n")
                R.utils::gzip(filename = temp.dist,
                    destname = paste0(outputFile, ".gz"), overwrite = TRUE)
            } else {
                data.table::fwrite(as.list(raw_lines), file = outputFile,
                                   sep = "\n")
            }
        }
        return(ret_obj)
    }

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")
    ret.df <- utils::read.table(text = ret.str[-1])
    ret.names <- ret.df[, 1]; ret.df <- ret.df[, -1]
    rownames(ret.df) <- ret.names; colnames(ret.df) <- ret.names

    if (!is.null(outputFile)) {
        if (compress) {
            temp.dist <- tempfile(fileext = ".dist"); on.exit(unlink(temp.dist), add = TRUE)
            data.table::fwrite(as.list(ret.str), file = temp.dist, sep = "\n")
            R.utils::gzip(filename = temp.dist,
                destname = paste0(outputFile, ".gz"), overwrite = TRUE)
        } else {data.table::fwrite(as.list(ret.str), file=outputFile, sep="\n")}
    }

    return(stats::as.dist(as.matrix(ret.df), diag = TRUE, upper = TRUE))
}

vcf2dist_parseWindowedLong <- function(lines) {
    if (length(lines) == 0L) {
        return(data.frame(chrom = character(), start = integer(),
                          end = integer(), sample_i = character(),
                          sample_j = character(), dist = numeric(),
                          stringsAsFactors = FALSE))
    }
    df <- utils::read.table(text = lines, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, comment.char = "")
    colnames(df)[1] <- "chrom"
    return(df)
}

vcf2dist_parseWindowedConcat <- function(lines) {
    out <- list()
    i <- 1L
    n <- length(lines)
    while (i <= n) {
        line <- lines[[i]]
        if (!startsWith(line, "# window")) {
            i <- i + 1L
            next
        }
        meta <- vcf2dist_parseWindowHeader(line)
        # next line: <nsamples>\t<nvariants>
        if (i + 1L > n) break
        dim_line <- lines[[i + 1L]]
        dims <- strsplit(dim_line, "\\s+")[[1]]
        nsamples <- as.integer(dims[1])
        # next nsamples lines: sample_name + nsamples distances
        body_start <- i + 2L
        body_end   <- i + 1L + nsamples
        if (body_end > n) break
        body <- stringr::str_replace_all(lines[body_start:body_end], "\t", " ")
        body_df <- utils::read.table(text = body, stringsAsFactors = FALSE)
        names_vec <- body_df[, 1]
        mat <- as.matrix(body_df[, -1])
        rownames(mat) <- names_vec
        colnames(mat) <- names_vec
        d <- stats::as.dist(mat, diag = TRUE, upper = TRUE)
        attr(d, "chrom") <- meta$chrom
        attr(d, "start") <- meta$start
        attr(d, "end") <- meta$end
        attr(d, "nvariants") <- meta$nvariants
        key <- paste0(meta$chrom, ":", meta$start, "-", meta$end)
        out[[key]] <- d
        i <- body_end + 1L
    }
    return(out)
}

vcf2dist_parseWindowHeader <- function(header) {
    chrom_m <- regmatches(header, regexpr("chrom=[^ ]+", header))
    start_m <- regmatches(header, regexpr("start=[0-9]+", header))
    end_m   <- regmatches(header, regexpr("end=[0-9]+", header))
    nvar_m  <- regmatches(header, regexpr("nvariants=[0-9]+", header))
    list(
        chrom = sub("chrom=", "", chrom_m),
        start = as.integer(sub("start=", "", start_m)),
        end   = as.integer(sub("end=", "", end_m)),
        nvariants = as.integer(sub("nvariants=", "", nvar_m))
    )
}

vcf2dist_checkParams <- function(inputFile, outputFile, threads, compress,
                                            verbose,
                                            windowBp, windowVariants,
                                            windowStep, windowMinVariants,
                                            longFormat) {
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

    if(!is.logical(compress)){
        stop("compress parameters must be logical.")
    }

    if (!is.numeric(threads) || (is.numeric(threads) && threads<1)) {
        stop("threads parameter must be positive integer.")
    }

    if (!is.logical(verbose)){
        stop("verbose",
             "must be logical.")
    }

    if (!is.null(windowBp) && !is.null(windowVariants)) {
        stop("windowBp and windowVariants are mutually exclusive.")
    }
    for (nm in c("windowBp", "windowVariants", "windowStep")) {
        v <- get(nm)
        if (!is.null(v) && (!is.numeric(v) || v < 1)) {
            stop(nm, " must be a positive integer.")
        }
    }
    if (!is.numeric(windowMinVariants) || windowMinVariants < 1) {
        stop("windowMinVariants must be a positive integer.")
    }
    if (!is.logical(longFormat)) {
        stop("longFormat must be logical.")
    }
    if (longFormat && is.null(windowBp) && is.null(windowVariants)) {
        stop("longFormat=TRUE requires windowBp or windowVariants to be set.")
    }
}
