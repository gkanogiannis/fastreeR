test_that("test parameters",{
    expect_error(vcf2clusters())
    expect_error(vcf2clusters(inputFile = 1))
    expect_error(vcf2clusters(inputFile = "thisdoesnotexist"))
    expect_error(vcf2clusters(inputFile = vcfFile, threads = 0))
    expect_error(vcf2clusters(inputFile = vcfFile, threads = "some"))
    expect_error(vcf2clusters(inputFile = vcfFile, ignoreMissing = 10))
    expect_error(vcf2clusters(inputFile = vcfFile, onlyHets = 1.1))
    expect_error(vcf2clusters(inputFile = vcfFile, ignoreHets = "some"))
    expect_error(vcf2clusters(inputFile = vcfFile, extra = 1))
    expect_error(vcf2clusters(inputFile = vcfFile, minClusterSize = 0))
    expect_error(vcf2clusters(inputFile = vcfFile, minClusterSize = "some"))
    expect_error(vcf2clusters(inputFile = vcfFile, cutHeight = -1))
    expect_error(vcf2clusters(inputFile = vcfFile, cutHeight = "some"))
})

test_that("test return is a list object with 3 items",{
    expect_type(vcf2clusters(inputFile = vcfFile), "list")
    expect_equal(length(vcf2clusters(inputFile = vcfFile)), 3)
})

test_that("test S1 is in the same cluster with S3",{
    expect_match(vcf2clusters(inputFile = vcfFile)[[3]][[2]], "S3.*S1|S1.*S3")
})

test_that("test with gzipped input",{
    expect_type(vcf2clusters(inputFile =
            system.file("extdata", "samples.vcf.gz", package = "fastreeR")),
                                                                        "list")
})
