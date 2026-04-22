#' Generate phylogenetic tree from samples of a VCF file
#'
#' This function calculates a distance matrix between the samples of a VCF file
#' as in \code{\link[fastreeR]{vcf2dist}}
#' and performs Hierarchical Clustering on this distance matrix
#' as in \code{\link[fastreeR]{dist2tree}}.
#' A phylogenetic tree is calculated by hierarchical clustering of the
#' distance matrix (complete linkage by default; single, complete, and
#' average linkage are supported by the Java backend).
#'
#' If the \code{bootstrap} parameter is set to a positive integer, the
#' Java backend performs streaming bootstrap sampling of variants for the
#' requested number of replicates. Bootstrap support values are encoded in
#' the returned Newick string at internal nodes (percent support across
#' replicates). Note that enabling bootstrapping increases runtime and
#' memory usage proportionally to the number of replicates.
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
#' \strong{Windowed mode.} Setting \code{windowBp} or \code{windowVariants}
#' (mutually exclusive) instructs the Java backend to emit one Newick tree per
#' genomic window. Windows never straddle chromosomes. Bootstrap is not
#' supported in windowed mode and will raise an error. The return value
#' becomes a \code{data.frame} with one row per window and columns
#' \code{chrom, start, end, nvariants, newick}.
#'
#' @param inputFile Input vcf file location (uncompressed or gzip compressed).
#' @param threads Number of java threads to use (default 1).
#' @param verbose Logical. If TRUE, enables verbose output from the Java backend.
#' @param bootstrap Number of bootstrap replicates to perform (default 0, no bootstrapping).
#' @param windowBp Optional positive integer. Emit one tree per window of N
#'   base pairs. Mutually exclusive with \code{windowVariants}.
#' @param windowVariants Optional positive integer. Emit one tree per N
#'   consecutive variants. Mutually exclusive with \code{windowBp}.
#' @param windowStep Optional positive integer. Window step (defaults to
#'   window size, i.e. tiled). Sliding windows are not yet implemented and
#'   will be rejected by the Java backend.
#' @param windowMinVariants Minimum number of variants required to emit a
#'   window (default 1; smaller windows are skipped silently).
#'
#' @return In non-windowed mode, a \code{\link[base]{character}} vector of the
#'   generated phylogenetic tree in Newick format. In windowed mode, a
#'   \code{data.frame} with columns \code{chrom, start, end, nvariants, newick}
#'   (one row per window).
#' @export
#'
#' @examples
#' my.tree <- vcf2tree(
#'     inputFile = system.file("extdata", "samples.vcf.gz",
#'         package = "fastreeR"
#'     )
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

vcf2tree <- function(inputFile, threads = 1,
                    verbose = FALSE, bootstrap = 0,
                    windowBp = NULL, windowVariants = NULL,
                    windowStep = NULL, windowMinVariants = 1L) {
    vcf2tree_checkParams(inputFile = inputFile, threads = threads,
                        verbose = verbose, bootstrap = bootstrap,
                        windowBp = windowBp, windowVariants = windowVariants,
                        windowStep = windowStep,
                        windowMinVariants = windowMinVariants)

    if (R.utils::isGzipped(inputFile)) {
        temp.in <- tempfile(fileext = ".vcf")
        on.exit(unlink(temp.in))
        R.utils::gunzip(filename = inputFile,
                        destname = temp.in,
                        remove = FALSE)
        inputFile <- temp.in
    }

    windowed <- !is.null(windowBp) || !is.null(windowVariants)

    bioinfojavautils <- rJava::.jnew(
        class="com/gkano/bioinfo/javautils/JavaUtils",
        class.loader = .rJava.class.loader
    )
    cmd_parts <- c(
        "VCF2TREE",
        "--numberOfThreads", threads,
        if (verbose) "--verbose" else "",
        "--input", inputFile,
        "--bootstrap", bootstrap
    )
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
    }
    cmd <- paste(cmd_parts, collapse = " ")

    temp.out <- tempfile(fileext = ".txt")
    on.exit(unlink(temp.out), add = TRUE)
    jSys <- rJava::J("java/lang/System")
    jOrigOut <- jSys$out
    jSys$setOut(rJava::.jnew("java/io/PrintStream", temp.out))
    bioinfojavautils$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    jSys$setOut(jOrigOut)

    if (windowed) {
        return(vcf2tree_parseWindowed(readLines(temp.out)))
    }

    ret.str <- stringr::str_replace_all(readLines(temp.out), "\t", " ")
    return(ret.str)
}

vcf2tree_parseWindowed <- function(lines) {
    chrom <- character(0); start <- integer(0); end <- integer(0)
    nvar  <- integer(0);   newick <- character(0)
    i <- 1L
    n <- length(lines)
    while (i <= n) {
        line <- lines[[i]]
        if (startsWith(line, "# window")) {
            meta <- vcf2dist_parseWindowHeader(line)
            tree_line <- ""
            j <- i + 1L
            while (j <= n && !startsWith(lines[[j]], "# window")) {
                if (nzchar(trimws(lines[[j]]))) {
                    tree_line <- lines[[j]]
                    break
                }
                j <- j + 1L
            }
            chrom  <- c(chrom,  meta$chrom)
            start  <- c(start,  meta$start)
            end    <- c(end,    meta$end)
            nvar   <- c(nvar,   meta$nvariants)
            newick <- c(newick, tree_line)
            i <- j + 1L
        } else {
            i <- i + 1L
        }
    }
    data.frame(chrom = chrom, start = start, end = end,
               nvariants = nvar, newick = newick,
               stringsAsFactors = FALSE)
}

vcf2tree_checkParams <- function(inputFile, threads, verbose, bootstrap,
                                 windowBp, windowVariants,
                                 windowStep, windowMinVariants) {
    if (!methods::is(inputFile, "character")){
        stop("inputFile must be a file location.")
    }

    if (is.null(inputFile) || !file.exists(inputFile)) {
        stop("inputFile=",inputFile," does not exist.")
    }

    if (!is.numeric(threads) || (is.numeric(threads) && threads<1)) {
        stop("threads parameter must be positive integer.")
    }

    if (!is.logical(verbose)){
        stop("verbose",
             "must be logical.")
    }

    if (!is.numeric(bootstrap) || (is.numeric(bootstrap) && bootstrap<0)) {
        stop("bootstrap parameter must be non-negative integer.")
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
    if (bootstrap > 0 && (!is.null(windowBp) || !is.null(windowVariants))) {
        stop("bootstrap is not supported in windowed mode (Java backend ",
             "rejects this combination).")
    }
}
