test_vcf2dist <- function() {
    RUnit::checkTrue(is.null(fastreeR::vcf2dist(inputfile = NULL)))
}
