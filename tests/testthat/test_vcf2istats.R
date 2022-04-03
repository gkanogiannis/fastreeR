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

test_that("test output file is correctly written with 4 lines",{
    istatsFile <- tempfile(fileext = ".istats")
    vcf2istats(inputFile = vcfFile, outputFile = istatsFile)
    expect_equal(length(readLines(istatsFile)), 4)
})
