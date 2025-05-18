
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="https://raw.githubusercontent.com/gkanogiannis/fastreeR/master/icon.png" alt="Project Icon" width="120"/>

# fastreeR: Fast Tree Reconstruction Tools for Genomics

<!-- badges: start -->

[![Bioconda](https://img.shields.io/conda/vn/bioconda/fastreer)](https://anaconda.org/bioconda/fastreer)
[![Docker
Pulls](https://img.shields.io/docker/pulls/gkanogiannis/fastreer)](https://hub.docker.com/r/gkanogiannis/fastreer)
[![PyPI
version](https://img.shields.io/pypi/v/fastreeR.svg)](https://pypi.org/project/fastreeR/)
<!-- badges: end -->

<!-- badges: start -->

BioC [![BioC
release](http://www.bioconductor.org/shields/build/release/bioc/fastreeR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/fastreeR)
<!-- badges: end -->

`fastreeR` is a hybrid toolkit combining a high-performance Java backend
([`BioInfoJava-Utils`](https://github.com/gkanogiannis/BioInfoJava-Utils)—a
modular Java library for bioinformatics pipelines) with flexible and
user-friendly interfaces across multiple platforms and environments,
enabling seamless integration into a variety of genomic workflows. It
enables fast computation of distance matrices and phylogenetic trees
from genetic variant data in **VCF** or genomic sequences in **FASTA**
format.

## Integration and Accessibility

`fastreeR` offers interface, which is accessible in the following ways:

- ✅ **Bioconda**: install with `conda install -c bioconda fastreer`
- ✅ **Docker**: available on
  [DockerHub](https://hub.docker.com/r/gkanogiannis/fastreer) and
  [GHCR](https://ghcr.io/gkanogiannis/fastreer) for containerized
  execution
- ✅ **PyPI**: install with `pip install fastreer`
- ✅ **Python CLI**: through a lightweight [Python
  wrapper](https://github.com/gkanogiannis/fastreeR/blob/devel/fastreeR.py)
  that calls the Java backend
- ✅ **R / Bioconductor**: via `rJava`
- ✅ **Galaxy**: Also available on Galaxy Toolshed.
- ✅ **Pure Java API**: developers can integrate this library directly
  in Java-based pipelines or software.

------------------------------------------------------------------------

- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation and Usage](#installation-and-usage)
- - [Conda](#via-conda)
  - [Docker](#via-docker)
  - [PyPI](#as-a-pypi-module)
  - [Python CLI](#via-a-python-cli-wrapper)
  - [R package](#as-an-r-package)
  - [Galaxy](#with-galaxy)
  - [From Java backend source](#from-java-backend-source)
- [Distances from VCF](#distances-from-vcf)
- [CLI Interface](#cli-interface)
  - [Commands](#commands)
  - [Examples](#examples)
  - [Options](#options-common-to-all-commands)
- [Integration with Java Backend](#integration-with-java-backend)
- [Integration with R](#integration-with-r)
- [Sample data](#sample-data)
- [Citation](#citation)
- [Author](#author)
- [License](#license)

------------------------------------------------------------------------

## Key Features

- ⚡ Ultra-fast computation of sample-wise cosine distances from large
  VCF and D2S k-mer based distances from FASTA files.
- 🌳 Generate agglomerative neighbor-joining phylogenetic trees directly
  from VCF or distance matrices.
- 🧵 Multithreaded execution for speed and scalability.
- Cluster distance matrices hierarchically with dynamic tree pruning.
- 🧰 Clean Python CLI for scripting and pipeline integration
- Streamlined integration with R via `rJava`
- Available on Galaxy Toolshed
- 🧬 Compatible with standard bioinformatics formats (PHYLIP, Newick)

------------------------------------------------------------------------

## Requirements

- Java 8+
- Python 3.7+
- Maven (if you want to build from the source)
- GNU/Linux, Windows or macOS

### Memory requirements for VCF input

At minimum, make sure to allocate for JVM at least 48 bytes per variant
per sample. If there are `n` samples and `m` variants allocate
`48 x n x m` bytes of RAM. For example, for processing a VCF file
containing data for 1 million variants and 1 thousand samples, allocate
at least : 48 x 10^6 x 10^3 = 48 x 10^9 bytes = 48GB of RAM. For optimal
execution, allocate more RAM than minimum. This will trigger less times
garbage collections and hence less pauses.

In order to allocate RAM, a special parameter needs to be passed while
JVM initializes. JVM parameters can be passed by setting
`java.parameters` option. The `-Xmx` parameter, followed (without space)
by an integer value and a letter, is used to tell JVM what is the
maximum amount of heap RAM that it can use. The letter in the parameter
(uppercase or lowercase), indicates RAM units. For example, parameters
`-Xmx1024m` or `-Xmx1024M` or `-Xmx1g` or `-Xmx1G`, all allocate 1
Gigabyte or 1024 Megabytes of maximum RAM for JVM.

In order to allocate 3GB of RAM for the JVM, through R code, use:

``` r
options(java.parameters = "-Xmx3G")
```

When using `fastreeR` as a CLI, then RAM allocation can be achieved with
the relevant argument `--mem MEM`.

A rough estimation for the required RAM, if sample and variant numbers
are not known, is half the size of the uncompressed VCF file. For
example for processing a VCF file, which uncompressed occupies 2GB of
disk space, allocate 1GB of RAM.

------------------------------------------------------------------------

## Installation and Usage

### Via Conda

``` bash
conda create -y -n fastreer-env -c bioconda fastreer && activate fastreer-env
fastreeR --help
```

### Via Docker

`fastreeR` is available as a lightweight, multithreaded,
platform-independent Docker image hosted on both **DockerHub** and
**GHCR**.

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

This: \* Mounts your working directory `$(pwd)` inside the container \*
Reads `input.vcf` and writes `output.nwk` relative to your host \* Uses
4 threads for faster computation

The Docker image includes: \* Java 17 \* Python3 \* All required `.jar`
libraries \* The `fastreeR.py` CLI entry point

Example: FASTA to distance

``` bash
docker run --rm -v $(pwd):/data gkanogiannis/fastreer \
    FASTA2DIST -i /data/sequences.fasta -o /data/sequences.dist -k 4 -t 2
```

Memory tuning Use the `--mem` option to control how much memory is
allocated to the Java backend:

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

This will install the fastreeR command-line tool (`fastreer`) and
include the Java backend jars required for running all commands.

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

Note: If you want to use a custom backend location, set the environment
variable `FASTREER_JAR_DIR`.

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

Search in Galaxy Tools for `fastreer` or ask your Galaxy Admin to
install it from toolshed.

### From java backend source

To build the Java backend from source code:

``` bash
git clone https://github.com/gkanogiannis/fastreeR.git
git clone https://github.com/gkanogiannis/BioInfoJava-Utils.git
pushd BioInfoJava-Utils
mvn clean initialize package && popd
```

Then copy the resulting `.jar` file(s) to the `fastreeR/inst/java/`
directory:

``` bash
cp BioInfoJava-Utils/bin/*.jar fastreeR/inst/java/
```

Finally run the tool from its Python CLI:

``` bash
python fastreeR/fastreeR.py
```

------------------------------------------------------------------------

## Distances from VCF

Calculates a cosine type dissimilarity measurement between the `n`
samples of a VCF file.

Biallelic or multiallelic (maximum 7 alternate alleles) SNP and/or INDEL
variants are considered, phased or not. Some VCF encoding examples are:

- heterozygous variants : `1/0` or `0/1` or `0/2` or `1|0` or `0|1` or
  `0|2`
- homozygous to the reference allele variants : `0/0` or `0|0`
- homozygous to the first alternate allele variants : `1/1` or `1|1`

If there are `n` samples and `m` variants, an `nxn` zero-diagonal
symmetric distance matrix is calculated. The calculated cosine type
distance (1-cosine_similarity)/2 is in the range `[0,1]` where value `0`
means completely identical samples (cosine is `1`), value `0.5` means
perpendicular samples (cosine is `0`) and value 1 means completely
opposite samples (cosine is `-1`).

The calculation is performed by a Java back-end implementation, that
supports multi-core CPU utilization and can be demanding in terms of
memory resources.

Output distances is a PHYLIP compatible file will contain `n+1` lines.
The first line contains the number `n` of samples and number `m` of
variants, separated by space. Each of the subsequent `n` lines contains
`n+1` values, separated by space. The first value of each line is a
sample name and the rest `n` values are the calculated distances of this
sample to all the samples. Example output file of the distances of 3
samples calculated from 1000 variants:

| 3 1000  |     |     |     |
|---------|-----|-----|-----|
| Sample1 | 0.0 | 0.5 | 0.2 |
| Sample2 | 0.5 | 0.0 | 0.9 |
| Sample3 | 0.2 | 0.9 | 0.0 |

------------------------------------------------------------------------

## CLI Interface

The Python CLI (`fastreeR.py`) interfaces with the Java backend via
`subprocess`, providing a unified command-line interface for all
supported tools.

### Commands

#### General Syntax

``` bash
python3 fastreeR.py <COMMAND> [OPTIONS]
```

| COMMAND      | Description                                      |
|--------------|--------------------------------------------------|
| `VCF2DIST`   | Compute a cosine distance matrix from a VCF file |
| `VCF2TREE`   | Compute a Newick NJ tree directly from a VCF     |
| `DIST2TREE`  | Compute a Newick NJ tree from a distance matrix  |
| `FASTA2DIST` | Compute a D2S distance matrix from a FASTA file  |

------------------------------------------------------------------------

### Examples

#### Compute Distance Matrix from VCF

``` bash
python fastreeR.py VCF2DIST -i input.vcf -o output.dist --threads 16 --verbose
```

#### Compute Newick NJ tree directly from a VCF file.

``` bash
python fastreeR.py VCF2TREE -i input.vcf -o output.nwk --threads 16 --verbose
```

#### Compute Tree from Distance Matrix

``` bash
python fastreeR.py DIST2TREE -i output.dist -o output.nwk
```

**Input format:** tab-separated PHYLIP-compatible matrix.

### Compute D2S k-mer distance matrix from a FASTA file.

``` bash
python3 fastreeR.py FASTA2DIST -i seqs.fasta -o output.dist -k 4 -t 2 --normalize
```

#### Pipe input from gzip-compressed file

``` bash
zcat input.vcf.gz | python fastreeR.py VCF2TREE -i - -o output.nwk
```

#### Print version and citation

``` bash
python fastreeR.py --version
```

### Output Examples

- Distance matrices: PHYLIP-compatible text
- Trees: Newick format
- Output is streamed line-by-line (suitable for large datasets)

------------------------------------------------------------------------

### Options (common to all commands)

- `-i, --input` : Input file (VCF or distance matrix). Use `-` for
  stdin.
- `-o, --output` : Output file. If omitted, prints to stdout.
- `-t, --threads` : Number of threads (default: 1).
- `--mem MEM` : Max RAM for JVM in GB (default: 1).
- `--lib LIB` : Path to the folder containing JAR libraries (default:
  inst/java)
- `--verbose` : Print progress information to stderr.
- `--pipe-stderr` : Pipe stderr and forward from Python (default: direct
  passthrough to terminal).
- `--version` : Print version and citation information.

------------------------------------------------------------------------

## Integration with Java Backend

The CLI wraps tools from the
[BioInfoJava-Utils](https://github.com/gkanogiannis/BioInfoJava-Utils)
project and dynamically builds the Java classpath from all `.jar` files
located in `inst/java/`.

------------------------------------------------------------------------

## Integration with R

All core functionality is available via the `fastreeR` R package
(Bioconductor/devel):

``` r
library(fastreeR)
tree <- vcf2tree("input.vcf")
plot(tree)
```

See [fastreeR R
manual](https://www.bioconductor.org/packages/release/bioc/manuals/fastreeR/man/fastreeR.pdf)
and [fastreeR R
vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/fastreeR/inst/doc/fastreeR_vignette.html)
for usage in R.

------------------------------------------------------------------------

## Sample data

Toy vcf, fasta and distance sample data files are provided in
`inst/extdata`.

### samples.vcf.gz

Sample VCF file of 100 individuals and 1000 variants, in Chromosome22,
from the 1K Genomes project. Original file available at
<http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/>

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

Sample FASTA file of 48 random bacteria RefSeq from
<ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/> .

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

> **Anestis Gkanogiannis (2016)**  
> *A scalable assembly-free variable selection algorithm for biomarker
> discovery from metagenomes*  
> BMC Bioinformatics 17, 311.  
> <https://doi.org/10.1186/s12859-016-1186-3>  
> <https://github.com/gkanogiannis/fastreeR>

------------------------------------------------------------------------

## Author

**Anestis Gkanogiannis**  
Website: <https://www.gkanogiannis.com>  
ORCID: [0000-0002-6441-0688](https://orcid.org/0000-0002-6441-0688)

------------------------------------------------------------------------

## License

`fastreeR` is licensed under the GNU General Public License v3.0.  
See the [LICENSE](LICENSE) file for details.

------------------------------------------------------------------------
