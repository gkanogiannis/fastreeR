# Perform Hierarchical Clustering and tree pruning on samples of VCF file

Performs Hierarchical Clustering on a distance matrix calculated as in
[`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
and generates a phylogenetic tree (complete linkage by default; single,
complete, and average linkage are supported by the Java backend), as in
[`dist2tree`](https://gkanogiannis.github.io/fastreeR/reference/dist2tree.md).
The phylogenetic tree is then pruned with
[`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)
to get clusters (as in
[`tree2clusters`](https://gkanogiannis.github.io/fastreeR/reference/tree2clusters.md)).

## Usage

``` r
vcf2clusters(
  inputFile,
  threads = 2,
  cutHeight = NULL,
  minClusterSize = 1,
  extra = TRUE,
  verbose = FALSE
)
```

## Arguments

- inputFile:

  Input vcf file location (uncompressed or gzip compressed).

- threads:

  Number of java threads to use.

- cutHeight:

  Define at which height to cut tree. Default automatically defined.

- minClusterSize:

  Minimum size of clusters. Default 1.

- extra:

  Boolean whether to use extra parameters for the
  [`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html).

- verbose:

  Logical. If TRUE, enables verbose output from the Java backend.

## Value

A list of :

- [`dist`](https://rdrr.io/r/stats/dist.html) distances object.

- [`character`](https://rdrr.io/r/base/character.html)` vector` of the
  generated phylogenetic tree in Newick format

- [`character`](https://rdrr.io/r/base/character.html)` vector` of the
  clusters. Each row contains data for a cluster, separated by space.
  The id of the cluster, the size of the cluster (number of elements)
  and the names of its elements, Cluster id 0 contains all the objects
  not assigned to a cluster (singletons). Example clusters output :

  |     |     |         |         |         |
  |-----|-----|---------|---------|---------|
  | 0   | 3   | Sample1 | Sample2 | Sample3 |
  | 1   | 3   | Sample4 | Sample5 | Sample6 |
  | 2   | 2   | Sample7 | Sample8 |         |
  | 3   | 2   | Sample9 | Sample0 |         |

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

The calculation is performed by a Java back-end implementation, that
supports multi-core CPU utilization and can be demanding in terms of
memory resources. By default a JVM is launched with a maximum memory
allocation of 512 MB. When this amount is not sufficient, the user needs
to reserve additional memory resources, before loading the package, by
updating the value of the `java.parameters` option. For example in order
to allocate 4GB of RAM, the user needs to issue
`options(java.parameters="-Xmx4g")` before
[`library(fastreeR)`](https://github.com/gkanogiannis/fastreeR).

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.clust <- vcf2clusters(
    inputFile = system.file("extdata", "samples.vcf.gz",
        package = "fastreeR"
    )
)
#>  ..cutHeight not given, setting it to 0.0407  ===>  99% of the (truncated) height range in dendro.
#>  ..done.
```
