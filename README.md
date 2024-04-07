
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastreeR

<!-- badges: start -->

[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/fastreeR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/fastreeR)
<!-- badges: end -->

The goal of fastreeR is to provide functions for calculating distance
matrix, building phylogenetic tree or performing hierarchical clustering
between samples, directly from a VCF or FASTA file.

## Requirements

A JDK, at least 8, is required and needs to be present before installing
`fastreeR`.

## Installation

To install `fastreeR` package:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("fastreeR")
```

You can install the development version of `fastreeR` like so:

``` r
devtools::install_github("gkanogiannis/fastreeR")
```

## Sample data

Toy vcf, fasta and distance sample data files are provided in
`inst/extdata`.

### samples.vcf.gz

Sample VCF file of 100 individuals and 1000 variants, in Chromosome22,
from the 1K Genomes project. Original file available at
<http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/>

``` r
vcfFile <- system.file("extdata", "samples.vcf.gz", package="fastreeR")
```

### samples.vcf.dist.gz

Distances from the previous sample VCF

``` r
vcfDist <- system.file("extdata", "samples.vcf.dist.gz", package="fastreeR")
```

### samples.vcf.istats

Individual statistics from the previous sample VCF

``` r
vcfIstats <- system.file("extdata", "samples.vcf.istats", package="fastreeR")
```

### samples.fasta.gz

Sample FASTA file of 48 random bacteria RefSeq from
<ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/> .

``` r
fastaFile <- system.file("extdata", "samples.fasta.gz", package="fastreeR")
```

### samples.fasta.dist.gz

Distances from the previous sample FASTA

``` r
fastaDist <- system.file("extdata", "samples.fasta.dist.gz",package="fastreeR")
```

## Memory requirements for VCF input

At minimum, make sure to allocate for JVM at least 10 bytes per variant
per sample. If there are `n` samples and `m` variants allocate
`10 x n x m` bytes of RAM. For example, for processing a VCF file
containing data for 1 million variants and 1 thousand samples, allocate
at least : 10^6 x 10^3 x 10 = 10^10 bytes = 10GB of RAM. For optimal
execution, allocate more RAM than minimum. This will trigger less times
garbage collections and hence less pauses.

In order to allocate RAM, a special parameter needs to be passed while
JVM initializes. JVM parameters can be passed by setting
`java.parameters` option. The `-Xmx` parameter, followed (without space)
by an integer value and a letter, is used to tell JVM what is the
maximum amount of heap RAM that it can use. The letter in the parameter
(uppercase or lowercase), indicates RAM units. For example, parameters
`-Xmx1024m` or `-Xmx1024M` or `-Xmx1g` or `-Xmx1G`, all allocate 1
Gigabyte or 1024 Megabytes of maximum RAM for JVM.

In order to allocate 3GB of RAM for the JVM, through R code, use:

``` r
options(java.parameters="-Xmx3G")
```

A rough estimation for the required RAM, if sample and variant numbers
are not known, is half the size of the uncompressed VCF file. For
example for processing a VCF file, which uncompressed occupies 2GB of
disk space, allocate 1GB of RAM.

## Distances from VCF

Calculates a cosine type dissimilarity measurement between the `n`
samples of a VCF file.

Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL
variants are considered, phased or not. Some VCF encoding examples are:

- heterozygous variants : `1/0` or `0/1` or `0/2` or `1|0` or `0|1` or
  `0|2`
- homozygous to the reference allele variants : `0/0` or `0|0`
- homozygous to the first alternate allele variants : `1/1` or `1|1`

If there are `n` samples and `m` variants, an `nxn` zero-diagonal
symmetric distance matrix is calculated. The calculated cosine type
distance (1-cosine_similarity)/2 is in the range `[0,1]` where value `0`
means completely identical samples (cosine is `1`), value `0.5` means
perpendicular samples (cosine is `0`) and value 1 means completely
opposite samples (cosine is `-1`).

The calculation is performed by a Java back-end implementation, that
supports multi-core CPU utilization and can be demanding in terms of
memory resources. By default a JVM is launched with a maximum memory
allocation of 512 MB. When this amount is not sufficient, the user needs
to reserve additional memory resources, before loading the package, by
updating the value of the `java.parameters` option. For example in order
to allocate 4GB of RAM, the user needs to issue
`options(java.parameters="-Xmx4g")` before `library(fastreeR)`.

Output file will contain `n+1` lines. The first line contains the number
`n` of samples and number `m` of variants, separated by space. Each of
the subsequent `n` lines contains `n+1` values, separated by space. The
first value of each line is a sample name and the rest `n` values are
the calculated distances of this sample to all the samples. Example
output file of the distances of 3 samples calculated from 1000 variants:

| 3 1000  |     |     |     |
|---------|-----|-----|-----|
| Sample1 | 0.0 | 0.5 | 0.2 |
| Sample2 | 0.5 | 0.0 | 0.9 |
| Sample3 | 0.2 | 0.9 | 0.0 |

## Distances from FASTA

## Tree from distances

## Clusters from tree
