# Perform Hierarchical Clustering and tree pruning on a phylogenetic tree

The phylogenetic tree is pruned with
[`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)
to get clusters.

## Usage

``` r
tree2clusters(
  treeStr,
  treeDistances = NULL,
  treeLabels = NULL,
  cutHeight = NULL,
  minClusterSize = 1,
  extra = TRUE,
  verbose = FALSE
)
```

## Arguments

- treeStr:

  A [`character`](https://rdrr.io/r/base/character.html)` vector` of a
  phylogenetic tree in Newick format

- treeDistances:

  `numeric `[`matrix`](https://rdrr.io/r/base/matrix.html) of distances,
  that were used to generate the tree. If NULL, it will be inferred from
  tree branch lengths.

- treeLabels:

  A [`character`](https://rdrr.io/r/base/character.html)` vector` of
  tree leaf labels.

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
my.clust <- tree2clusters(
    treeStr = dist2tree(
        inputDist = system.file("extdata", "samples.vcf.dist.gz",
            package = "fastreeR"
        )
    ),
    verbose = TRUE
)
#>  ..cutHeight not given, setting it to 0.0793  ===>  99% of the (truncated) height range in dendro.
#>  ..done.
```
