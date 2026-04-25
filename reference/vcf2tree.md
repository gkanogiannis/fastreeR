# Generate phylogenetic tree from samples of a VCF file

This function calculates a distance matrix between the samples of a VCF
file as in
[`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
and performs Hierarchical Clustering on this distance matrix as in
[`dist2tree`](https://gkanogiannis.github.io/fastreeR/reference/dist2tree.md).
A phylogenetic tree is calculated by hierarchical clustering of the
distance matrix (complete linkage by default; single, complete, and
average linkage are supported by the Java backend).

## Usage

``` r
vcf2tree(
  inputFile,
  threads = 1,
  verbose = FALSE,
  bootstrap = 0,
  windowBp = NULL,
  windowVariants = NULL,
  windowStep = NULL,
  windowMinVariants = 1L
)
```

## Arguments

- inputFile:

  Input vcf file location (uncompressed or gzip compressed).

- threads:

  Number of java threads to use (default 1).

- verbose:

  Logical. If TRUE, enables verbose output from the Java backend.

- bootstrap:

  Number of bootstrap replicates to perform (default 0, no
  bootstrapping).

- windowBp:

  Optional positive integer. Emit one tree per window of N base pairs.
  Mutually exclusive with `windowVariants`.

- windowVariants:

  Optional positive integer. Emit one tree per N consecutive variants.
  Mutually exclusive with `windowBp`.

- windowStep:

  Optional positive integer. Window step (defaults to window size, i.e.
  tiled). Sliding windows are not yet implemented and will be rejected
  by the Java backend.

- windowMinVariants:

  Minimum number of variants required to emit a window (default 1;
  smaller windows are skipped silently).

## Value

In non-windowed mode, a
[`character`](https://rdrr.io/r/base/character.html) vector of the
generated phylogenetic tree in Newick format. In windowed mode, a
`data.frame` with columns `chrom, start, end, nvariants, newick` (one
row per window).

## Details

If the `bootstrap` parameter is set to a positive integer, the Java
backend performs streaming bootstrap sampling of variants for the
requested number of replicates. Bootstrap support values are encoded in
the returned Newick string at internal nodes (percent support across
replicates). Note that enabling bootstrapping increases runtime and
memory usage proportionally to the number of replicates.

Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL
variants are considered, phased or not. Some VCF encoding examples are:

- heterozygous variants : `1/0` or `0/1` or `0/2` or `1|0` or `0|1` or
  `0|2`

- homozygous to the reference allele variants : `0/0` or `0|0`

- homozygous to the first alternate allele variants : `1/1` or `1|1`

If there are `n` samples and `m` variants, an `nxn` zero-diagonal
symmetric distance matrix is calculated. The calculated cosine type
distance (1-cosine_similarity)/2 is in the range \[0,1\] where value 0
means completely identical samples (cosine is 1), value 0.5 means
perpendicular samples (cosine is 0) and value 1 means completely
opposite samples (cosine is -1).

The calculation is performed by a Java backend implementation, that
supports multi-core CPU utilization and can be demanding in terms of
memory resources. By default a JVM is launched with a maximum memory
allocation of 512 MB. When this amount is not sufficient, the user needs
to reserve additional memory resources, before loading the package, by
updating the value of the `java.parameters` option. For example in order
to allocate 4GB of RAM, the user needs to issue
`options(java.parameters="-Xmx4g")` before
[`library(fastreeR)`](https://github.com/gkanogiannis/fastreeR).

**Windowed mode.** Setting `windowBp` or `windowVariants` (mutually
exclusive) instructs the Java backend to emit one Newick tree per
genomic window. Windows never straddle chromosomes. Bootstrap is not
supported in windowed mode and will raise an error. The return value
becomes a `data.frame` with one row per window and columns
`chrom, start, end, nvariants, newick`.

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.tree <- vcf2tree(
    inputFile = system.file("extdata", "samples.vcf.gz",
        package = "fastreeR"
    )
)
```
