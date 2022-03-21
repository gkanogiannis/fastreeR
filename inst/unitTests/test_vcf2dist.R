test_vcf2dist <- function() {
    checkTrue(is.na(fastreeR::vcf2dist(inputfile = NULL)))
}
