.onLoad <- function(libname, pkgname) {
    # Check that Java version is at least 8
    rJava::.jinit()
    jv <- rJava::.jcall("java/lang/System", "S",
                        "getProperty", "java.runtime.version")
    if(substr(jv, 1L, 2L) == "1.") {
        jvn <- as.numeric(
                    paste0(
                        strsplit(jv, "[.]")[[1L]][seq_len(2)], collapse = "."
                    )
                )
        if(jvn < 1.8)
            stop("Java >= 8 is needed for this package but not available")
    }

    rJava::.jpackage(
        name = "fastreeR",
        own.loader = TRUE,
        parameters = c(getOption("java.parameters"),
                        "-XX:+UseG1GC",
                        "-XX:+UseStringDeduplication")
    )
}
