test_that("test parameters",{
    expect_error(vcf2dist())
    expect_error(vcf2dist(inputFile = 1))
    expect_error(vcf2dist(inputFile = "thisdoesnotexist"))
    expect_error(vcf2dist(inputFile = vcfFile, outputFile = 1))
    expect_error(vcf2dist(inputFile = vcfFile, outputFile = ""))
    expect_error(vcf2dist(inputFile = vcfFile, threads = 0))
    expect_error(vcf2dist(inputFile = vcfFile, threads = "some"))
    expect_error(vcf2dist(inputFile = vcfFile, ignoreMissing = 10))
    expect_error(vcf2dist(inputFile = vcfFile, onlyHets = 1.1))
    expect_error(vcf2dist(inputFile = vcfFile, ignoreHets = "some"))
    expect_error(vcf2dist(inputFile = vcfFile, compress = 1))
})

test_that("test return is a dist object",{
    expect_s3_class(vcf2dist(inputFile = vcfFile), "dist")
})

test_that("test S1 is closer (less distance) to S3 that to S2",{
    my.dist <- vcf2dist(inputFile = vcfFile)
    expect_true(my.dist[2]<my.dist[1])
})

test_that("test S1 identical (distance=0) to S3 when ignoreHets=TRUE",{
    expect_equal(vcf2dist(inputFile = vcfFile, ignoreHets = TRUE)[2], 0)
})

test_that("test S2 identical (distance=0) to S3 when onlyHets=TRUE",{
    expect_equal(vcf2dist(inputFile = vcfFile, onlyHets = TRUE)[3], 0)
})
