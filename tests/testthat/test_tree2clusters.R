test_that("test parameters",{
    expect_error(tree2clusters())
    expect_error(tree2clusters(treeStr = 1))
    expect_error(tree2clusters(treeStr = treeStr, extra = 1))
    expect_error(tree2clusters(treeStr = treeStr, minClusterSize = 0))
    expect_error(tree2clusters(treeStr = treeStr, minClusterSize = "some"))
    expect_error(tree2clusters(treeStr = treeStr, cutHeight = -1))
    expect_error(tree2clusters(treeStr = treeStr, cutHeight = "some"))
    expect_error(tree2clusters(treeStr = treeStr, treeDistances = 1))
    expect_error(tree2clusters(treeStr = treeStr, treeDistances = "some"))
    expect_error(tree2clusters(treeStr = treeStr, treeLabels = 1))
})

test_that("test return is a character vector",{
    expect_type(tree2clusters(treeStr = treeStr), "character")
})

test_that("test S1 is in the same cluster with S3",{
    expect_match(tree2clusters(treeStr = treeStr)[[2]], "S3.*S1|S1.*S3")
})
