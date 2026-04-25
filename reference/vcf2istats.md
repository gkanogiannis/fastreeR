# Calculate various per sample statistics from a VCF file

Only biallelic SNPs are considered. For each sample, the following
statistics are calculated :

- INDIV : Sample name

- N_SITES : Total number of SNPs

- N_HET : Number of SNPs with heterozygous call (`0/1` or `0|1` or `1/0`
  or `1|0`)

- N_ALT : Number of SNPs with alternate homozygous call (`1/1` or `1|1`)

- N_REF : Number of SNPs with reference homozygous call (`0/0` or `0|0`)

- N_MISS : Number of SNPs with missing call (`./.` or `.|.`)

- P_HET : Percentage of heterozygous calls

- P_ALT : Percentage of alternate homozygous calls

- P_REF : Percentage of reference homozygous calls

- P_MISS : Percentage of missing calls (missing rate)

## Usage

``` r
vcf2istats(inputFile, outputFile = NULL)
```

## Arguments

- inputFile:

  Input vcf file location (uncompressed or gzip compressed).

- outputFile:

  Output samples statistics file location.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) of sample
statistics.

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.istats <- vcf2istats(
    inputFile =
        system.file("extdata", "samples.vcf.gz", package = "fastreeR")
)
```
