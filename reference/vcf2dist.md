# Calculate distances between samples of a VCF file

This function calculates a cosine type dissimilarity measurement between
the `n` samples of a VCF file.

## Usage

``` r
vcf2dist(
  inputFile,
  outputFile = NULL,
  threads = 2,
  compress = FALSE,
  verbose = FALSE,
  windowBp = NULL,
  windowVariants = NULL,
  windowStep = NULL,
  windowMinVariants = 1L,
  longFormat = FALSE
)
```

## Arguments

- inputFile:

  Input vcf file location (uncompressed or gzip compressed).

- outputFile:

  Output distances file location.

- threads:

  Number of java threads to use.

- compress:

  Compress output (adds .gz extension).

- verbose:

  Logical. If TRUE, enables verbose output from the Java backend.

- windowBp:

  Optional positive integer. Emit one matrix per window of N base pairs.
  Mutually exclusive with `windowVariants`.

- windowVariants:

  Optional positive integer. Emit one matrix per N consecutive variants.
  Mutually exclusive with `windowBp`.

- windowStep:

  Optional positive integer. Window step (defaults to window size, i.e.
  tiled). Sliding windows are not yet implemented and will be rejected
  by the Java backend.

- windowMinVariants:

  Minimum number of variants required to emit a window (default 1;
  smaller windows are skipped silently).

- longFormat:

  Logical. In windowed mode, return a long-form `data.frame` instead of
  a list of `dist` objects.

## Value

In non-windowed mode, a [`dist`](https://rdrr.io/r/stats/dist.html)
distances object. In windowed mode, either a named `list` of `dist`
objects (default) or a long-form `data.frame` (when
`longFormat = TRUE`).

## Details

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

Output file, if provided, will contain `n+1` lines. The first line
contains the number `n` of samples and number `m` of variants, separated
by space. Each of the subsequent `n` lines contains `n+1` values,
separated by space. The first value of each line is a sample name and
the rest `n` values are the calculated distances of this sample to all
the samples. Example output file of the distances of 3 samples
calculated from 1000 variants:

|        |     |         |     |
|--------|-----|---------|-----|
| 3 1000 |     | Sample1 | 0.0 |
| 0.5    | 0.2 | Sample2 | 0.5 |
| 0.0    | 0.9 | Sample3 | 0.2 |

**Windowed mode.** Setting `windowBp` or `windowVariants` (mutually
exclusive) instructs the Java backend to emit one distance matrix per
genomic window. Windows never straddle chromosomes. The return value
changes accordingly:

- `longFormat = FALSE` (default): a named `list` of
  [`dist`](https://rdrr.io/r/stats/dist.html) objects, one per window.
  List names follow the format `"chrom:start-end"`.

- `longFormat = TRUE`: a single `data.frame` with columns
  `chrom, start, end, sample_i, sample_j, dist`.

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.dist <- vcf2dist(
    inputFile = system.file("extdata", "samples.vcf.gz",
        package = "fastreeR"
    )
)
```
