# Perform Hierarchical Clustering and tree pruning on a distance matrix

Performs Hierarchical Clustering on a distance matrix (i.e. calculated
with
[`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
or
[`fasta2dist`](https://gkanogiannis.github.io/fastreeR/reference/fasta2dist.md))
and generates a phylogenetic tree (complete linkage by default; single,
complete, and average linkage are supported by the Java backend), as in
[`dist2tree`](https://gkanogiannis.github.io/fastreeR/reference/dist2tree.md).
The phylogenetic tree is then pruned with
[`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)
to get clusters (as in
[`tree2clusters`](https://gkanogiannis.github.io/fastreeR/reference/tree2clusters.md)).

## Usage

``` r
dist2clusters(
  inputDist,
  cutHeight = NULL,
  minClusterSize = 1,
  extra = TRUE,
  verbose = FALSE
)
```

## Arguments

- inputDist:

  Input distances file location (generated with
  [`vcf2dist`](https://gkanogiannis.github.io/fastreeR/reference/vcf2dist.md)
  or
  [`fasta2dist`](https://gkanogiannis.github.io/fastreeR/reference/fasta2dist.md)).
  File can be gzip compressed. Or a
  [`dist`](https://rdrr.io/r/stats/dist.html) distances object.

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

## References

Java implementation: <https://github.com/gkanogiannis/BioInfoJava-Utils>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
my.clust <- dist2clusters(
    inputDist =
        system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR"),
    verbose = TRUE
)
#>  ..cutHeight not given, setting it to 0.0793  ===>  99% of the (truncated) height range in dendro.
#>  ..done.
```
