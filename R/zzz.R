.onLoad <- function(libname, pkgname) {
    java_params <- getOption("java.parameters")
    if (is.null(java_params)) java_params <- character(0)

    rJava::.jpackage(
        name = pkgname,
        lib.loc = libname,
        own.loader = TRUE,
        parameters = c(
            "-Djava.awt.headless=true",
            "-XX:+UseG1GC",
            "-XX:+UseStringDeduplication",
            "-XX:+UseCompressedOops",
            "-Dfile.encoding=UTF-8",
            java_params
        )
    )
}
