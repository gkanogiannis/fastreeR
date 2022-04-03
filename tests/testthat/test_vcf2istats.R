test_that("test parameters",{
    expect_error(vcf2istats())
    expect_error(vcf2istats(inputFile = 1))
    expect_error(vcf2istats(inputFile = "thisdoesnotexist"))
    expect_error(vcf2istats(inputFile = vcfFile, outputFile = 1))
    expect_error(vcf2istats(inputFile = vcfFile, outputFile = ""))
})

test_that("test return is a data.frame object with 10 columns",{
    expect_s3_class(vcf2istats(inputFile = vcfFile), "data.frame")
    expect_equal(length(vcf2istats(inputFile = vcfFile)), 10)
})
