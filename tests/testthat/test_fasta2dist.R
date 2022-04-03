test_that("test parameters",{
    expect_error(fasta2dist())
    expect_error(fasta2dist(1))
    expect_error(fasta2dist("thisdoesnotexist"))
    expect_error(fasta2dist(fastaFile, outputFile = 1))
    expect_error(fasta2dist(fastaFile, outputFile = ""))
    expect_error(fasta2dist(fastaFile, threads = 0))
    expect_error(fasta2dist(fastaFile, threads = "some"))
    expect_error(fasta2dist(fastaFile, kmer = 0))
    expect_error(fasta2dist(fastaFile, kmer = "some"))
    expect_error(fasta2dist(fastaFile, normalize = 10))
    expect_error(fasta2dist(fastaFile, normalize = "some"))
    expect_error(fasta2dist(fastaFile, compress = 10))
    expect_error(fasta2dist(fastaFile, compress = "some"))
})

test_that("test return is a dist object",{
    expect_s3_class(fasta2dist(fastaFile), "dist")
})

test_that("test S1 is closer (less distance) to S2 that to S3",{
    my.dist <- fasta2dist(fastaFile)
    expect_true(my.dist[1]<my.dist[2])
})
