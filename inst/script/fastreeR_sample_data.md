## Sample data

Toy vcf, fasta and distance sample data files are provided in `inst/extdata`.

### samples.vcf.gz
Sample VCF file of 100 individuals and 1000 variants, in Chromosome22, 
from the 1K Genomes project. Original file available at
[http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/)

```
vcfFile <- system.file("extdata", "samples.vcf.gz", package="fastreeR")
```

### samples.vcf.dist.gz
Distances from the previous sample VCF
```
vcfDist <- system.file("extdata", "samples.vcf.dist.gz", package="fastreeR")
```

### samples.vcf.istats
Individual statistics from the previous sample VCF
```
vcfIstats <- system.file("extdata", "samples.vcf.istats", package="fastreeR")
```

### samples.fasta.gz
Sample FASTA file of 48 random bacteria RefSeq from 
[ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/) .
```
fastaFile <- system.file("extdata", "samples.fasta.gz", package="fastreeR")
```

### samples.fasta.dist.gz
Distances from the previous sample FASTA
```
fastaDist <- system.file("extdata", "samples.fasta.dist.gz",package="fastreeR")
```
