
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastreeR

<!-- badges: start -->
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

Toy sample files are provided in `inst/extdata`.

Sample VCF file of 100 samples and 1000 variants, in Chromosome22, from
the 1K Genomes project. Original file available at
<http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/>

``` r
vcf.file <- system.file("extdata", "samples.vcf.gz", package="fastreeR")
```

Distances from sample VCF

``` r
vcf.dist <- system.file("extdata", "samples.vcf.dist.gz", package="fastreeR")
```

Sample FASTA file of 48 random bacteria RefSeq from
<ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/> .

``` r
fasta.file <- system.file("extdata", "samples.fasta.gz", package="fastreeR")
```

Distances from sample FASTA

``` r
fasta.dist <- system.file("extdata", "samples.fasta.dist.gz",package="fastreeR")
```

## Memory requirements for VCF input

At minimum, make sure to allocate for JVM at least 10 bytes per variant
per sample. If there are samples and variants allocate 10xnxm bytes of
RAM. For example, for processing a VCF file containing data for 1
million variants and 1 thousand samples, allocate at least : 10^6 x 10^3
x 10 = 10^10 bytes = 10GB of RAM. For optimal execution, allocate more
RAM than minimum. This will trigger less times garbage collections and
hence less pauses.

In order to allocate 3GB of RAM for the JVM, through R code, use:

``` r
options(java.parameters="-Xmx3G")
```

A rough estimation for the required RAM, if sample and variant numbers
are not known, is half the size of the uncompressed VCF file. For
example for processing a VCF file, which uncompressed occupies 2GB of
disk space, allocate 1GB of RAM.

## Distances from VCF

Calculates a cosine type dissimilarity measurement between the samples
of a VCF file.

Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL
variants are considered, phased or not. Some VCF encoding examples are:

If there are samples and variants, an zero-diagonal symmetric distance
matrix is calculated. The calculated cosine type distance
(1-cosine_similarity)/2 is in the range \[0,1\] where value 0 means
completely identical samples (cosine is 1), value 0.5 means
perpendicular samples (cosine is 0) and value 1 means completely
opposite samples (cosine is -1).

The calculation is performed by a Java backend implementation, that
supports multi-core CPU utilization and can be demanding in terms of
memory resources. By default a JVM is launched with a maximum memory
allocation of 512 MB. When this amount is not sufficient, the user needs
to reserve additional memory resources, before loading the package, by
updating the value of the option. For example in order to allocate 4GB
of RAM, the user needs to issue before .

Output file will contain lines. The first line contains the number of
samples and number of variants, separated by space. Each of the
subsequent lines contains values, separated by space. The first value of
each line is a sample name and the rest values are the calculated
distances of this sample to all the samples. Example output file of the
distances of 3 samples calculated from 1000 variants:

## Distances from FASTA

## Tree from distances

## Clusters from tree
