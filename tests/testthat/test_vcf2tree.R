test_that("test parameters",{
    expect_error(vcf2tree())
    expect_error(vcf2tree(inputFile = 1))
    expect_error(vcf2tree(inputFile = "thisdoesnotexist"))
    expect_error(vcf2tree(inputFile = vcfFile, threads = 0))
    expect_error(vcf2tree(inputFile = vcfFile, threads = "some"))
    expect_error(vcf2tree(inputFile = vcfFile, ignoreMissing = 10))
    expect_error(vcf2tree(inputFile = vcfFile, onlyHets = 1.1))
    expect_error(vcf2tree(inputFile = vcfFile, ignoreHets = "some"))
})

test_that("test return is a character vector",{
    expect_type(vcf2tree(inputFile = vcfFile), "character")
})


test_that("test with gzipped input",{
    expect_type(vcf2tree(inputFile =
            system.file("extdata", "samples.vcf.gz", package = "fastreeR")),
                                                                    "character")
})
