test_that("test parameters",{
    expect_error(dist2tree())
    expect_error(dist2tree(inputDist = 1))
    expect_error(dist2tree(inputDist = "thisdoesnotexist"))
})

test_that("test return is a character vector",{
    expect_type(dist2tree(inputDist = distFile), "character")
})

test_that("test with gzipped input",{
    expect_type(dist2tree(inputDist =
        system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")),
                                                                    "character")

})

test_that("test with dist input",{
    expect_type(dist2tree(inputDist = vcf2dist(inputFile = vcfFile)),
                                                                    "character")
})
