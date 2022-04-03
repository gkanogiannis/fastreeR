#Dummy vcf (3 samples, 4 variants) temp file
#S1 is identical to S3 when ignoring heterozygous calls
#S2 is identical to S3 when using only heterozygous calls
#S1 is closer (less distance) to S3 that to S2

vcfFile <- tempfile(fileext = ".vcf")
vcfStr <-paste0("##fileformat=VCFv4.1\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n",
        "1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\n",
        "1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\n",
        "1\t3\t.\tT\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\n",
        "1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1\n")
write(vcfStr, file = vcfFile, sep = "")


#Dummy dist temp file (generated from the dummy vcf file)

distFile <- tempfile(fileext = ".dist")
vcf2dist(inputFile = vcfFile, outputFile = distFile, compress = FALSE)


#Dummy tree in newick format

treeStr <- vcf2tree(inputFile = vcfFile)


#Dummy fasta (3 samples of 10bp length) temp file

fastaFile <- tempfile(fileext = ".fasta")
fastaStr <- paste0(
    ">S1\n",
    "ACGTACGTAA\n",
    ">S2\n",
    "ACGTACGTCC\n",
    ">S3\n",
    "AAAAAGGGGG\n"
)
write(fastaStr, file = fastaFile, sep = "")
