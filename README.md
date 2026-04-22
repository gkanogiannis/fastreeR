
<img src="https://raw.githubusercontent.com/gkanogiannis/fastreeR/master/icon.png" alt="Project Icon" width="120"/>

# fastreeR: Fast Tree Reconstruction Tools for Genomics (VCF/FASTA to Distance/Tree)

<!-- badges: start -->
![Conda Version](https://img.shields.io/conda/v/bioconda/fastreer)![Conda Downloads](https://img.shields.io/conda/dn/bioconda/fastreer) ![Docker Version](https://img.shields.io/docker/v/gkanogiannis/fastreer?label=docker)![Docker Pulls](https://img.shields.io/docker/pulls/gkanogiannis/fastreer?label=pulls) ![PyPI Version](https://img.shields.io/pypi/v/fastreer)![PyPI Downloads](https://img.shields.io/pypi/dm/fastreer)
![Galaxy Version](https://img.shields.io/galaxytoolshed/v/fastreer_vcf2tree/gkanogiannis/fastreer_vcf2tree)![Galaxy Downloads](https://img.shields.io/galaxytoolshed/downloads/fastreer_vcf2tree/gkanogiannis) [![Bioc history](https://bioconductor.org/shields/years-in-bioc/fastreeR.svg)](https://bioconductor.org/packages/release/bioc/html/fastreeR.html#since)[![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/fastreeR.svg)](http://bioconductor.org/packages/stats/bioc/fastreeR/)
<!-- badges: end -->

`fastreeR` is a hybrid toolkit combining a high-performance Java backend ([`BioInfoJava-Utils`](https://github.com/gkanogiannis/BioInfoJava-Utils), a modular Java library for bioinformatics pipelines) with flexible and user-friendly interfaces across multiple platforms and environments, enabling seamless integration into a variety of genomic workflows.
It enables fast computation of distance matrices and phylogenetic trees from genetic variant data in **VCF** or genomic sequences in **FASTA** format.

## Integration and Accessibility

`fastreeR` offers interface, which is accessible in the following ways:

* 🆕 **Java Backend ([v2.7.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.7.0)) !!** introduces **windowed / streaming VCF distance & tree output**. Emit one distance matrix (or Newick tree) per genomic window of N base pairs (`--window-bp`) or per N consecutive variants (`--window-variants`) for `VCF2DIST` and `VCF2TREE`, with optional long-form TSV output (`--long`). Windows never straddle chromosomes.
* Java Backend ([v2.5.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.5.0)) introduces **embedding-based distance calculation** for VCF files. Provide pre-computed variant embeddings (from genomic language models like [BioFM](https://huggingface.co/m42-health/BioFM-265M), DNA-BERT, Nucleotide Transformer, etc.) to weight variant contributions during distance computation.
* Java Backend ([v2.3.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.3.0)) supports reading from gzip (for example .gz), bzip2 (for example .bz2) and xz compressed VCF files.
* Java Backend ([v2.2.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.2.0)) implements streaming bootstrap; from VCF file get a newick tree with encoded bootstrap support values.
* Java Backend ([v2.0.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/2.0.0)) 100x times **FAST**re**ER** and only a couple hundred MB RAM needed. Java 11+ suggested.
* **Bioconda**: install with `conda install -c bioconda fastreer` ([recipe](https://bioconda.github.io/recipes/fastreer/README.html))
* **Docker**: available on [DockerHub](https://hub.docker.com/r/gkanogiannis/fastreer) and [GHCR](https://ghcr.io/gkanogiannis/fastreer) for containerized execution
* **PyPI**: install with `pip install fastreer` ([repository](https://pypi.org/project/fastreer/))
* **Python CLI**: through a lightweight [Python wrapper](https://github.com/gkanogiannis/fastreeR/blob/devel/fastreeR.py) that calls the Java backend
* **R / Bioconductor**: via `rJava` ([package](https://bioconductor.org/packages/fastreeR/))
* **Galaxy**: available on Galaxy [Toolshed](https://toolshed.g2.bx.psu.edu/view/gkanogiannis/fastreer/26013530719e).
* **Pure Java API**: developers can integrate this library directly in Java-based pipelines or software.

------------------------------------------------------------------------

- [fastreeR: Fast Tree Reconstruction Tools for Genomics (VCF/FASTA to Distance/Tree)](#fastreer-fast-tree-reconstruction-tools-for-genomics-vcffasta-to-distancetree)
  - [Integration and Accessibility](#integration-and-accessibility)
  - [Key Features](#key-features)
  - [Requirements](#requirements)
    - [Memory requirements for VCF input](#memory-requirements-for-vcf-input)
  - [Installation and Usage](#installation-and-usage)
    - [Via Conda](#via-conda)
    - [Via Docker](#via-docker)
    - [As a PyPI Module](#as-a-pypi-module)
    - [Via a Python CLI wrapper](#via-a-python-cli-wrapper)
    - [As an R package](#as-an-r-package)
    - [With Galaxy](#with-galaxy)
    - [From java backend source](#from-java-backend-source)
  - [Distances from VCF](#distances-from-vcf)
  - [Embedding-Based Distance Calculation](#embedding-based-distance-calculation)
    - [How It Works](#how-it-works)
    - [Embedding File Formats](#embedding-file-formats)
    - [Embedding Command Line Options](#embedding-command-line-options)
    - [Embedding Examples](#embedding-examples)
  - [Windowed / Streaming Output](#windowed--streaming-output)
    - [How Windowing Works](#how-windowing-works)
    - [Windowing Command Line Options](#windowing-command-line-options)
    - [Output Formats](#output-formats)
    - [Windowing Examples](#windowing-examples)
    - [Windowing Limitations](#windowing-limitations)
    - [Windowed output from R](#windowed-output-from-r)
  - [CLI Interface](#cli-interface)
    - [Commands](#commands)
      - [General Syntax](#general-syntax)
    - [Examples](#examples)
      - [Compute Distance Matrix from VCF](#compute-distance-matrix-from-vcf)
      - [Compute Newick tree directly from a VCF file.](#compute-newick-tree-directly-from-a-vcf-file)
      - [Compute Tree from Distance Matrix](#compute-tree-from-distance-matrix)
      - [Compute D2S k-mer distance matrix from a FASTA file.](#compute-d2s-k-mer-distance-matrix-from-a-fasta-file)
      - [Generate Variant Embeddings from VCF using BioFM](#generate-variant-embeddings-from-vcf-using-biofm)
      - [Pipe input from gzip-compressed file](#pipe-input-from-gzip-compressed-file)
    - [Output Examples](#output-examples)
    - [Options (common to all commands)](#options-common-to-all-commands)
    - [Embedding options (VCF2DIST and VCF2TREE only)](#embedding-options-vcf2dist-and-vcf2tree-only)
    - [Windowing options (VCF2DIST and VCF2TREE only)](#windowing-options-vcf2dist-and-vcf2tree-only)
    - [VCF2EMB options (embedding generation)](#vcf2emb-options-embedding-generation)
  - [Integration with Java Backend](#integration-with-java-backend)
  - [Integration with R](#integration-with-r)
  - [Sample data](#sample-data)
    - [samples.vcf.gz](#samplesvcfgz)
    - [samples.vcf.dist.gz](#samplesvcfdistgz)
    - [samples.vcf.istats](#samplesvcfistats)
    - [samples.fasta.gz](#samplesfastagz)
    - [samples.fasta.dist.gz](#samplesfastadistgz)
  - [Citation](#citation)
  - [Author](#author)
  - [License](#license)
  
------------------------------------------------------------------------

## Key Features

* 📁 Input from standard VCF (gz, bzip2, xz compressed or uncompressed) and FASTA files.
* 🪟 **Windowed / streaming output** emits one distance matrix or Newick tree per genomic window (by base pairs or variant count) for `VCF2DIST` and `VCF2TREE`.
* 🧠 **Embedding-based distance calculation** using pre-computed variant embeddings from genomic language models.
* 🥾 Streaming bootstrap support from VCF to NEWICK.
* 🚀 With a superior multithreaded concurrency model and minimal RAM usage, from GBs down to just MBs!
* ⚡ Ultra-fast computation of sample-wise cosine distances from large VCF and D2S k-mer based distances from FASTA files.
* Generate phylogenetic trees directly from VCF or distance matrices using **hierarchical clustering** (single, complete, or average linkage; complete by default).
* Multithreaded execution for speed and scalability.
* Cluster distance matrices hierarchically with dynamic tree pruning.
* Clean Python CLI for scripting and pipeline integration
* Streamlined integration with R via `rJava`
* Available on Galaxy Toolshed
* Compatible with standard bioinformatics formats (PHYLIP, Newick)

------------------------------------------------------------------------

## Requirements

* Java 11+
* Python 3.7+
* Maven (if you want to build from the source)
* GNU/Linux, Windows or macOS

### Memory requirements for VCF input

**No more GBs of RAM!** Only the distance matrix is kept in memory:

* `4 bytes x (#samples²) x #threads`
* Example: 1000 samples with 32 threads → **\~128MB RAM**

**VCF caching is minimal:**
Only **2 VCF lines per thread** are pre-cached.

* In the simple diploid case (e.g., `0/1`, `1|0`), each genotype requires \~4 characters (8 bytes).
* For 1000 samples and 32 threads, this adds up to **\~1MB RAM**.

JVM will need at least 64-128 MB in order to efficiently run.

**Total memory footprint: just a few hundred MB, even for large datasets.**

~~It is not straightforward to define a strict minimum amount of RAM required for a given number of SNPs and samples, as JVM behavior can vary across different systems and configurations.
From our own experiments, a rough estimate for the minimum usable memory is around 10 bytes per variant per sample.
For example, a VCF file with 1 million variants and 1,000 samples would require at least 10 x 10⁶ x 10³ = 10 GB of allocated memory.
However, running with this minimal allocation may result in frequent and prolonged garbage collection events, leading to significantly longer runtimes.
For optimal execution, we recommend allocating 15-20 bytes per variant per sample (i.e., 15-20 GB for the same example), which reduces garbage collection overhead and ensures smoother performance.~~

In order to allocate RAM, a special parameter needs to be passed while JVM initializes. JVM parameters can be passed by setting `java.parameters` option.
The `-Xmx` parameter, followed (without space) by an integer value and a letter, is used to tell JVM what is the maximum amount of heap RAM that it can use.
The letter in the parameter (uppercase or lowercase), indicates RAM units.
For example, parameters `-Xmx1024m` or `-Xmx1024M` or `-Xmx1g` or `-Xmx1G`, allocate 1 Gigabyte or 1024 Megabytes of maximum RAM for JVM.

In order to allocate 1024MB of RAM for the JVM, through R code, use:

``` r
options(java.parameters = "-Xmx1024M")
```

When using `fastreeR` as a CLI, then RAM allocation in MB can be achieved with the relevant argument `--mem MEM`.

------------------------------------------------------------------------

## Installation and Usage

### Via Conda

`fastreeR` is available on Bioconda. You can install it in a new conda environment like so:

``` bash
conda create -y -n fastreer-env -c bioconda fastreer && activate fastreer-env
fastreeR --help
```

### Via Docker

`fastreeR` is available as a lightweight, multithreaded, platform-independent Docker image hosted on both **DockerHub** and **GHCR**.

From DockerHub:

``` bash
docker pull gkanogiannis/fastreer:latest
```

Or from GitHub Container Registry (GHCR):

``` bash
docker pull ghcr.io/gkanogiannis/fastreer:latest
```

To compute a tree directly from a VCF file:

``` bash
docker run --rm -v $(pwd):/data gkanogiannis/fastreer:latest \
    VCF2TREE -i /data/input.vcf -o /data/output.nwk --threads 4
```

This:

* Mounts your working directory `$(pwd)` inside the container
* Reads `input.vcf` and writes `output.nwk` relative to your host
* Uses 4 threads for faster computation

The Docker image includes:

* Java 21
* Python3
* All required `.jar` libraries
* The `fastreeR.py` CLI entry point

Example: FASTA to distance

``` bash
docker run --rm -v $(pwd):/data gkanogiannis/fastreer \
    FASTA2DIST -i /data/sequences.fasta -o /data/sequences.dist -k 4 -t 2
```

Memory tuning.
Use the `--mem` option to control how much memory is allocated to the Java backend:

``` bash
docker run --rm -v $(pwd):/data gkanogiannis/fastreer \
    VCF2TREE -i /data/input.vcf -o /data/output.nwk --mem 128
```

> Internally, this sets the Java heap to `-Xmx128G`.

### As a PyPI Module

You can install the Python CLI directly from PyPI using:

``` bash
pip install fastreer
```

This will install the fastreeR command-line tool (`fastreer`) and include the Java backend jars required for running all commands.

To check it installed correctly:

``` bash
fastreeR --version
```

### Via a Python CLI wrapper

Another easy method for using `fastreeR` is by its Python CLI:

``` bash
git clone https://github.com/gkanogiannis/fastreeR.git
python fastreeR/fastreeR.py
```

Note: If you want to use a custom backend location, set the environment variable `FASTREER_JAR_DIR`.

### As an R package

To install `fastreeR` as an R package:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("fastreeR")
```

You can install the development version of `fastreeR` R package like so:

``` r
devtools::install_github("gkanogiannis/fastreeR")
```

### With Galaxy

Search in Galaxy Tools for `fastreer` or ask your Galaxy Admin to install it from toolshed.

### From java backend source

To build the Java backend from source code:

``` bash
git clone https://github.com/gkanogiannis/fastreeR.git
git clone https://github.com/gkanogiannis/BioInfoJava-Utils.git
pushd BioInfoJava-Utils
mvn clean initialize package && popd
```

Then copy the resulting `.jar` file(s) to the `fastreeR/inst/java/` directory:

``` bash
cp BioInfoJava-Utils/bin/*.jar fastreeR/inst/java/
```

Finally run the tool from its Python CLI:

``` bash
python fastreeR/fastreeR.py
```

------------------------------------------------------------------------

## Distances from VCF

Calculates a cosine type dissimilarity measurement between the `n` samples of a VCF file.

Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL variants are considered, phased or not. Some VCF encoding examples are:

* heterozygous variants : `1/0` or `0/1` or `0/2` or `1|0` or `0|1` or `0|2`
* homozygous to the reference allele variants : `0/0` or `0|0`
* homozygous to the first alternate allele variants : `1/1` or `1|1`

If there are `n` samples and `m` variants, an `nxn` zero-diagonal symmetric distance matrix is calculated.
The calculated cosine type distance (1-cosine_similarity)/2 is in the range `[0,1]` where value `0` means completely identical samples (cosine is `1`), value `0.5` means perpendicular samples (cosine is `0`) and value 1 means completely opposite samples (cosine is `-1`).

The calculation is performed by a Java back-end implementation, that supports multi-core CPU utilization and can be demanding in terms of memory resources.

Output distances is a PHYLIP compatible file will contain `n+1` lines.
The first line contains the number `n` of samples and number `m` of variants, separated by space.
Each of the subsequent `n` lines contains `n+1` values, separated by space.
The first value of each line is a sample name and the rest `n` values are the calculated distances of this sample to all the samples.
Example output file of the distances of 3 samples calculated from 1000 variants:

| 3 1000  |     |     |     |
|---------|-----|-----|-----|
| Sample1 | 0.0 | 0.5 | 0.2 |
| Sample2 | 0.5 | 0.0 | 0.9 |
| Sample3 | 0.2 | 0.9 | 0.0 |

------------------------------------------------------------------------

## Embedding-Based Distance Calculation

Version 2.5.0 of the Java backend introduces support for **embedding-based distance calculation** in `VCF2DIST` and `VCF2TREE`. This feature allows you to incorporate pre-computed variant embeddings (e.g., from genomic language models like [BioFM](https://huggingface.co/m42-health/BioFM-265M), DNA-BERT, Nucleotide Transformer, or custom embeddings) to compute distances in embedding space rather than genotype space.

### How It Works

Instead of computing cosine similarity directly from genotype vectors, the embedding mode:

1. Projects each sample into embedding space: `H_i = Σ_v dosage_i^v × e_v`
2. Computes cosine distance between sample embeddings

This captures functional relationships between variants - samples with alleles at functionally similar positions become more similar in embedding space.

### Embedding File Formats

**TSV Format:**

```tsv
#VARIANT_ID  DIM_0   DIM_1   DIM_2   ...
chr1:12345:A:G  0.123   -0.456  0.789   ...
chr1:67890:C:T  0.567   0.123   -0.890  ...
```

**HuggingFace JSON Format:**

```json
{
  "model_name": "genomic-model-name",
  "embedding_dim": 768,
  "variants": [
    {"id": "chr1:12345:A:G", "embedding": [0.123, -0.456, ...]},
    {"id": "chr1:67890:C:T", "embedding": [0.567, 0.123, ...]}
  ]
}
```

### Embedding Command Line Options

| Option                 | Description                                                                     |
|------------------------|---------------------------------------------------------------------------------|
| `-e, --embeddings`     | Path to variant embeddings file                                                 |
| `--embeddings-format`  | Format: `TSV` or `HUGGINGFACE` (auto-detected if not specified)                 |
| `--variant-key`        | Variant key format: `CHROM_POS`, `CHROM_POS_REF_ALT` (default), or `VCF_ID`     |

### Embedding Examples

``` bash
# Distance matrix with embeddings (TSV format, auto-detected)
python fastreeR.py VCF2DIST -i samples.vcf.gz -o distances.tsv -e variant_embeddings.tsv -t 4

# Tree with embeddings and bootstrap (HuggingFace format)
python fastreeR.py VCF2TREE -i samples.vcf.gz -o tree.nwk -e embeddings.json --embeddings-format HUGGINGFACE -b 100

# Standard mode (no embeddings) - existing behavior
python fastreeR.py VCF2DIST -i samples.vcf.gz -o distances.tsv
```

Variants without matching embeddings are automatically skipped, and the tool reports how many variants were used vs. skipped.

------------------------------------------------------------------------

## Windowed / Streaming Output

Version 2.7.0 of the Java backend introduces **windowed output** for `VCF2DIST` and `VCF2TREE`. Instead of producing a single genome-wide distance matrix or tree, the tools can stream one matrix (or Newick tree) per genomic window. This enables local-ancestry analyses, introgression scans, recombination-rate studies, and any workflow that needs sample relationships measured along the genome.

### How Windowing Works

Variants are streamed in input order and grouped into windows defined either by base-pair span (`--window-bp`) or by consecutive variant count (`--window-variants`). When a window closes, all worker threads synchronize on a barrier, the per-window distance matrix is reduced from shared accumulators, the writer emits the window, and the accumulators are zeroed before the next window opens. **Windows never straddle chromosomes**; a contig change always closes the current window.

The non-windowed code path is unchanged and remains byte-identical to previous releases.

### Windowing Command Line Options

| Option              | Description                                                                                                |
|---------------------|------------------------------------------------------------------------------------------------------------|
| `--window-bp`       | Emit one matrix/tree per window of N base pairs (mutually exclusive with `--window-variants`)              |
| `--window-variants` | Emit one matrix/tree per N consecutive variants (mutually exclusive with `--window-bp`)                    |
| `--step`            | Window step. Defaults to window size (tiled). Sliding windows (`step != size`) are not yet implemented.    |
| `--min-variants`    | Minimum number of variants required to emit a window (default 1; smaller windows are skipped silently)     |
| `--long`            | (`VCF2DIST` only) Emit long-form TSV `chrom, start, end, sample_i, sample_j, dist` instead of matrices     |

### Output Formats

`VCF2DIST` default (concatenated matrices), one block per window:

```text
# window chrom=chr1 start=0 end=100000 nvariants=842 nsamples=3
3	842
s1	0	0.4231	0.5102
s2	0.4231	0	0.3987
s3	0.5102	0.3987	0
# window chrom=chr1 start=100000 end=200000 nvariants=917 nsamples=3
...
```

`VCF2DIST --long`, single TSV with one row per sample pair per window:

```text
chrom	start	end	sample_i	sample_j	dist
chr1	0	100000	s1	s2	0.4231
chr1	0	100000	s1	s3	0.5102
chr1	0	100000	s2	s3	0.3987
...
```

`VCF2TREE`, one Newick tree per window, prefixed by a header comment:

```text
# window chrom=chr1 start=0 end=100000 nvariants=842 nsamples=3
(s1:0.21,(s2:0.19,s3:0.18):0.05);
# window chrom=chr1 start=100000 end=200000 nvariants=917 nsamples=3
(s2:0.20,(s1:0.22,s3:0.17):0.04);
...
```

### Windowing Examples

``` bash
# Distance matrices in 100kb tiled windows
python fastreeR.py VCF2DIST -i samples.vcf.gz -o per_window.dist --window-bp 100000 -t 4

# Long-form TSV, one matrix per 500 consecutive variants
python fastreeR.py VCF2DIST -i samples.vcf.gz -o per_window.tsv --window-variants 500 --long -t 4

# Per-window phylogenetic trees (Newick)
python fastreeR.py VCF2TREE -i samples.vcf.gz -o per_window.nwk --window-bp 250000 -t 4

# Skip windows with fewer than 50 variants
python fastreeR.py VCF2DIST -i samples.vcf.gz -o per_window.dist --window-bp 100000 --min-variants 50
```

### Windowing Limitations

- **Sliding windows** (`--step` different from window size) are reserved for a future release; passing them throws an error.
- **Bootstrap** (`-b` / `--bootstrap`) is rejected when combined with windowing.
- **Embeddings** (`-e` / `--embeddings`) are rejected when combined with windowing.

### Windowed output from R

`vcf2dist()` and `vcf2tree()` accept the same windowing parameters
(`windowBp`, `windowVariants`, `windowStep`, `windowMinVariants`, plus
`longFormat` for `vcf2dist`). When any window parameter is set the return
value changes to one of:

* `vcf2dist(..., windowBp = 100000)` — named `list` of `dist` objects, one per window (names are `"chrom:start-end"`).
* `vcf2dist(..., windowVariants = 500, longFormat = TRUE)` — single long-form `data.frame` with columns `chrom, start, end, sample_i, sample_j, dist`.
* `vcf2tree(..., windowBp = 250000)` — `data.frame` with columns `chrom, start, end, nvariants, newick`.

``` r
library(fastreeR)
vcf <- system.file("extdata", "samples.vcf.gz", package = "fastreeR")

# Per-window distance matrices (list of dist)
windows <- vcf2dist(vcf, windowBp = 100000)
length(windows); head(names(windows))

# Per-window trees as a data.frame
trees <- vcf2tree(vcf, windowVariants = 500)
trees[1, ]
```

------------------------------------------------------------------------

## CLI Interface

The Python CLI (`fastreeR.py`) interfaces with the Java backend via `subprocess`, providing a unified command-line interface for all supported tools.

### Commands

#### General Syntax

``` bash
python3 fastreeR.py <COMMAND> [OPTIONS]
```

| COMMAND      | Description                                                      |
|--------------|------------------------------------------------------------------|
| `VCF2DIST`   | Compute a cosine distance matrix from a VCF file (genome-wide or per window) |
| `VCF2TREE`   | Compute a Newick hierarchical-clustering tree from a VCF (genome-wide or per window) |
| `DIST2TREE`  | Compute a Newick hierarchical-clustering tree from a distance matrix |
| `FASTA2DIST` | Compute a D2S distance matrix from a FASTA file                  |
| `VCF2EMB`    | Generate variant embeddings from VCF using BioFM language model  |

------------------------------------------------------------------------

### Examples

#### Compute Distance Matrix from VCF

``` bash
python fastreeR.py VCF2DIST -i input.vcf -o output.dist --threads 16 --verbose
```

#### Compute Newick tree directly from a VCF file.

``` bash
python fastreeR.py VCF2TREE -i input.vcf -o output.nwk --threads 16 --verbose
```

You can also request bootstrap replicates directly from the VCF source. The Java backend will perform streaming bootstrap sampling and encode bootstrap support values at internal nodes of the returned Newick string.
For example:

``` bash
python fastreeR.py VCF2TREE -i input.vcf -o output_with_boot.nwk --threads 8 --bootstrap 100
```

The generated Newick will contain node support values (percentage across replicates) which can be inspected with phylogenetic tools such as `ape` in R.

#### Compute Tree from Distance Matrix

``` bash
python fastreeR.py DIST2TREE -i output.dist -o output.nwk
```

**Input format:** tab-separated PHYLIP-compatible matrix.

#### Compute D2S k-mer distance matrix from a FASTA file.

``` bash
python3 fastreeR.py FASTA2DIST -i seqs.fasta -o output.dist -k 4 -t 2 --normalize
```

#### Generate Variant Embeddings from VCF using BioFM

The `VCF2EMB` command uses the [BioFM-265M](https://huggingface.co/m42-health/BioFM-265M) genomic language model to generate embeddings for each variant in a VCF file. These embeddings can then be used with `VCF2DIST` or `VCF2TREE` for embedding-based distance calculation.

**Supports gzipped input files:** VCF (`.vcf.gz`), reference genome (`.fa.gz`, `.fasta.gz`, `.fna.gz`), and annotation (`.gff.gz`, `.gff3.gz`) files are automatically decompressed during processing.

**Prerequisites:**

1. Python 3.11 environment (required by biofm-eval):
   ```bash
   conda create -n fastreer-env python=3.11
   conda activate fastreer-env
   ```

2. Install PyTorch:
   ```bash
   pip install torch  # CPU only
   # Or with CUDA: pip install torch --index-url https://download.pytorch.org/whl/cu121
   ```

3. Install biofm-eval from source (not available on PyPI):
   ```bash
   git clone https://github.com/m42-health/biofm-eval.git
   cd biofm-eval
   pip install -e .
   ```

4. Download reference genome (GRCh38): [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)

5. Download gene annotations (GENCODE v38): [GENCODE](https://www.gencodegenes.org/human/release_38.html)

``` bash
# Generate embeddings in TSV format (supports gzipped inputs)
python fastreeR.py VCF2EMB -i input.vcf.gz -o embeddings.tsv \
    -r GRCh38.fna.gz -a gencode.v38.annotation.gff3.gz --verbose

# Generate embeddings in HuggingFace JSON format
python fastreeR.py VCF2EMB -i input.vcf.gz -o embeddings.json \
    -r GRCh38.fna -a gencode.v38.annotation.gff3 -f HUGGINGFACE

# Use GPU for faster processing
python fastreeR.py VCF2EMB -i input.vcf.gz -o embeddings.tsv \
    -r GRCh38.fna -a gencode.v38.annotation.gff3 --device cuda

# Process only first 1000 variants
python fastreeR.py VCF2EMB -i input.vcf.gz -o embeddings.tsv \
    -r GRCh38.fna -a gencode.v38.annotation.gff3 --max-variants 1000
```

You can set default paths via environment variables:
``` bash
export BIOFM_REFERENCE_GENOME=/path/to/GRCh38.fna.gz
export BIOFM_GENE_ANNOTATION=/path/to/gencode.v38.annotation.gff3.gz
python fastreeR.py VCF2EMB -i input.vcf.gz -o embeddings.tsv
```

#### Pipe input from gzip-compressed file

``` bash
zcat input.vcf.gz | python fastreeR.py VCF2TREE -i - -o output.nwk
```

### Output Examples

* Distance matrices: PHYLIP-compatible text
* Trees: Newick format
* Output is streamed line-by-line (suitable for large datasets)

------------------------------------------------------------------------

### Options (common to all commands)

* `-i, --input` : Input file (VCF or distance matrix). Use `-` for stdin.
* `-o, --output` : Output file. If omitted, prints to stdout.
* `-t, --threads` : Number of threads (default: 1).
* `--mem MEM` : Max RAM for JVM in MB (default: 256).
* `--lib LIB` : Path to the folder containing backend JAR libraries (default: inst/java)
* `--verbose` : Print progress information to stderr.
* `--pipe-stderr` : Pipe stderr and forward from Python (default: direct passthrough to terminal).
* `--version` : Print version and citation information.

### Embedding options (VCF2DIST and VCF2TREE only)

* `-e, --embeddings` : Path to variant embeddings file for embedding-based distance calculation.
* `--embeddings-format` : Embeddings file format: `TSV` or `HUGGINGFACE` (auto-detected if not specified).
* `--variant-key` : Variant key format for embedding lookup: `CHROM_POS`, `CHROM_POS_REF_ALT` (default), or `VCF_ID`.

### Windowing options (VCF2DIST and VCF2TREE only)

* `--window-bp N` : Emit one matrix/tree per window of `N` base pairs (mutually exclusive with `--window-variants`).
* `--window-variants N` : Emit one matrix/tree per `N` consecutive variants (mutually exclusive with `--window-bp`).
* `--step N` : Window step (defaults to window size, i.e. tiled). Sliding windows are not yet implemented.
* `--min-variants N` : Minimum number of variants required to emit a window (default 1).
* `--long` : (`VCF2DIST` only) Emit long-form TSV `chrom, start, end, sample_i, sample_j, dist` instead of concatenated matrices.

### VCF2EMB options (embedding generation)

* `-i, --input` : Input VCF file.
* `-o, --output` : Output embeddings file (default: stdout).
* `-r, --reference` : Path to reference genome FASTA file (or set `BIOFM_REFERENCE_GENOME` env var).
* `-a, --annotation` : Path to gene annotation GFF3 file (or set `BIOFM_GENE_ANNOTATION` env var).
* `-m, --model` : HuggingFace model name or local path (default: `m42-health/BioFM-265M`).
* `-f, --format` : Output format: `TSV` or `HUGGINGFACE` (default: `TSV`).
* `--variant-key` : Variant key format in output: `CHROM_POS`, `CHROM_POS_REF_ALT` (default), or `VCF_ID`.
* `--max-variants` : Maximum number of variants to process (default: all).
* `--batch-size` : Batch size for embedding extraction (default: 32).
* `--device` : Device for model inference: `cuda` or `cpu` (default: auto-detect).
 
------------------------------------------------------------------------

## Integration with Java Backend

The CLI wraps tools from the [BioInfoJava-Utils](https://github.com/gkanogiannis/BioInfoJava-Utils) project and dynamically builds the Java classpath from all `.jar` files located in `inst/java/`.

------------------------------------------------------------------------

## Integration with R

All core functionality is available via the `fastreeR` R package (Bioconductor/devel):

``` r
library(fastreeR)
tree <- vcf2tree("input.vcf")
plot(tree)
```

See [fastreeR R manual](https://www.bioconductor.org/packages/release/bioc/manuals/fastreeR/man/fastreeR.pdf) and [fastreeR R vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/fastreeR/inst/doc/fastreeR_vignette.html) for usage in R.

------------------------------------------------------------------------

## Sample data

Toy vcf, fasta and distance sample data files are provided in `inst/extdata`.

### samples.vcf.gz

Sample VCF file of 100 individuals and 1000 variants, in Chromosome22, from the 1K Genomes project. Original file available at [http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/)

``` r
vcfFile <- system.file("extdata", "samples.vcf.gz", package = "fastreeR")
```

### samples.vcf.dist.gz

Distances from the previous sample VCF

``` r
vcfDist <- system.file("extdata", "samples.vcf.dist.gz", package = "fastreeR")
```

### samples.vcf.istats

Individual statistics from the previous sample VCF

``` r
vcfIstats <- system.file("extdata", "samples.vcf.istats", package = "fastreeR")
```

### samples.fasta.gz

Sample FASTA file of 48 random bacteria RefSeq from [ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/).

``` r
fastaFile <- system.file("extdata", "samples.fasta.gz", package = "fastreeR")
```

### samples.fasta.dist.gz

Distances from the previous sample FASTA

``` r
fastaDist <- system.file("extdata", "samples.fasta.dist.gz", package = "fastreeR")
```

------------------------------------------------------------------------

## Citation

If you use `fastreeR` in your research, please cite:

> Anestis Gkanogiannis (2016)
> _A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes_  
> _BMC Bioinformatics_ 17, 311.  
> <https://doi.org/10.1186/s12859-016-1186-3>  
> <https://github.com/gkanogiannis/fastreeR>

------------------------------------------------------------------------

## Author

Anestis Gkanogiannis  
Bioinformatics/ML Scientist  
Linkedin: [https://www.linkedin.com/in/anestis-gkanogiannis/](https://www.linkedin.com/in/anestis-gkanogiannis/)  
Website: [https://github.com/gkanogiannis](https://github.com/gkanogiannis)  
ORCID: [0000-0002-6441-0688](https://orcid.org/0000-0002-6441-0688)

------------------------------------------------------------------------

## License

`fastreeR` is licensed under the GNU General Public License v3.0.  
See the [LICENSE](LICENSE.md) file for details.

------------------------------------------------------------------------
