# Generate phylogenetic tree from samples of a distance matrix

Performs Hierarchical Clustering on a distance matrix (i.e. calculated
with
[`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
or
[`fasta2dist`](https://gkanogiannis.github.io/fastreeR/reference/fasta2dist.md))
and generates a phylogenetic tree (complete linkage by default; single,
complete, and average linkage are supported by the Java backend).

## Usage

``` r
dist2tree(inputDist, verbose = FALSE)
```

## Arguments

- inputDist:

  Input distances file location (generated with
  [`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
  or
  [`fasta2dist`](https://gkanogiannis.github.io/fastreeR/reference/fasta2dist.md)).
  File can be gzip compressed. Or a
  [`dist`](https://rdrr.io/r/stats/dist.html) distances object.

- verbose:

  Logical. If TRUE, enables verbose output from the Java backend.

## Value

A [`character`](https://rdrr.io/r/base/character.html)` vector` of the
generated phylogenetic tree in Newick format.

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.tree <- dist2tree(
    inputDist =
    system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR"),
    verbose = TRUE
)
```
