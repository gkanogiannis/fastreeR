test_that("vcf2tree returns Newick with bootstrap values when requested", {
    skip_on_cran()
    vcf <- system.file("extdata", "samples.vcf.gz", package = "fastreeR")
    # use a small number of replicates to keep test fast
    tr_str <- vcf2tree(inputFile = vcf, threads = 1, bootstrap = 10)
    expect_type(tr_str, "character")
    # basic Newick check: must end with semicolon
    expect_true(grepl(";\\s*$", tr_str))
    # bootstrap values are typically numbers in internal node labels, e.g. :0.95 or [95]
    # Check for presence of numbers between parentheses/brackets or node labels
    expect_true(grepl("[0-9]{1,3}", tr_str))
})
