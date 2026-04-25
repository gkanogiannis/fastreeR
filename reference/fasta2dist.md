# Calculate distances between sequences of a FASTA file

This function calculates a d2_S type dissimilarity measurement between
the `n` sequences (which can represent samples) of a FASTA file. See
[doi:10.1186/s12859-016-1186-3](https://doi.org/10.1186/s12859-016-1186-3)
for more details.

## Usage

``` r
fasta2dist(
  ...,
  outputFile = NULL,
  threads = 2,
  kmer = 6,
  normalize = FALSE,
  compress = TRUE,
  verbose = FALSE
)
```

## Arguments

- ...:

  Input fasta files locations (uncompressed or gzip compressed).

- outputFile:

  Output distances file location.

- threads:

  Number of java threads to use.

- kmer:

  Kmer length to use for analyzing fasta sequences.

- normalize:

  Normalize on sequences length.

- compress:

  Compress output (adds .gz extension).

- verbose:

  Logical. If TRUE, enables verbose output from the Java backend.

## Value

A [`dist`](https://rdrr.io/r/stats/dist.html) distances object of the
calculation.

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.dist <- fasta2dist(
    inputfile = system.file("extdata", "samples.fasta.gz",
        package = "fastreeR"
    )
)
```
