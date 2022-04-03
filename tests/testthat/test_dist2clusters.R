test_that("test parameters",{
    expect_error(dist2clusters())
    expect_error(dist2clusters(inputDist = 1))
    expect_error(dist2clusters(inputDist = "thisdoesnotexist"))
    expect_error(dist2clusters(inputDist = distFile, extra = 1))
    expect_error(dist2clusters(inputDist = distFile, minClusterSize = 0))
    expect_error(dist2clusters(inputDist = distFile, minClusterSize = "some"))
    expect_error(dist2clusters(inputDist = distFile, cutHeight = -1))
    expect_error(dist2clusters(inputDist = distFile, cutHeight = "some"))
})

test_that("test return is a list object with 3 items",{
    expect_type(dist2clusters(inputDist = distFile), "list")
    expect_equal(length(dist2clusters(inputDist = distFile)), 2)
})

test_that("test S1 is in the same cluster with S3",{
    expect_match(dist2clusters(inputDist = distFile)[[2]][[2]], "S3.*S1|S1.*S3")
})
