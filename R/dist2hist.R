#' Plot a histogram of distance matrix
#'
#' Plots a histogram of distances, from a distance matrix
#' (i.e. calculated with \code{\link[fastreeR]{vcf2dist}}).
#'
#' @param input.dist Input distances file location
#' (generated with \code{\link[fastreeR]{vcf2dist}}).
#' File can be gzip compressed.
#' Or a \code{\link[stats]{dist}} distances object.
#' @param outputfile Output histogram png image file location.
#' @param settings Settings for histogram. bins,width,height
#'
#' @usage
#' dist2hist(
#'     input.dist,
#'     outputfile = NULL,
#'     settings = "100,1024,768"
#' )
#'
#' @return A \code{\link[png]{readPNG} nativeRaster} of the produced histogram.
#' @export
#'
#' @examples
#' my.hist <- dist2hist(
#'     input.dist =
#'         system.file("extdata", "samples.dist.gz", package = "fastreeR")
#' )
#' grid::grid.raster(my.hist)
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references Java implementation:
#' \url{https://github.com/gkanogiannis/BioInfoJava-Utils}

dist2hist <- function(
                    input.dist,
                    outputfile = NULL,
                    settings = "100,1024,768") {
    if (is.null(input.dist) ||
        (!methods::is(input.dist, "dist") &&
            !methods::is(input.dist, "character")) ||
        (methods::is(input.dist, "character") && !file.exists(input.dist))) {
        return(NA)
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
        utils::write.table(
            as.matrix(input.dist), file = temp.in, sep = " ", row.names = TRUE,
            col.names = FALSE, append = TRUE
        )
        inputfile <- temp.in
    }

    bioinfojavautils <- rJava::J(
        class="ciat/agrobio/javautils/JavaUtils",
        class.loader = .rJava.class.loader)
    temp.out <- tempfile(fileext = ".png")
    on.exit(unlink(temp.out))
    cmd <- paste(
        "DIST2Hist",
        inputfile, "--output", temp.out,
        "--settings", settings, sep = " "
    )

    bioinfojavautils$main(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))
    gc()
    rJava::J("java.lang.Runtime")$getRuntime()$gc()

    if (!is.null(outputfile)) {
        return(png::readPNG(temp.out, native = TRUE))
    } else {
        return(NA)
    }
}
