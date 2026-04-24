# fastreeR NEWS

## fastreeR 2.2.0 with java backend 2.7.1

- Windowed / streaming VCF distance and tree output. Emit one distance matrix (or Newick tree) per genomic window of N base pairs or per N consecutive variants for VCF2DIST and VCF2TREE.

## fastreeR 2.1.3 with java backend 2.5.0 (2026-02-02)

- Added embedding-based distance calculation for VCF files (backend 2.5.0).

## fastreeR 2.1.2 with java backend 2.3.0 (2025-12-01)

- Implement reading compressed VCF input (gz, bzip2, xz).

## fastreeR 2.1.0 with java backend 2.2.0 (2025-11-01)

- Major update: enhanced streaming, faster bootstrapping, and backend v2.2.0 integration.

## Version 2.0.0 (2025-10-29)

- Implements streaming bootstrap; from VCF to newick with support values.
- Streaming ultra-fast VCF processing; almost zero RAM needed.
- Update backend (2.2.0) and CLI. Preparing for major version upgrade.

## Version 1.99.0 (2025-05-29)

- Update backend and CLI. Preparing for major version upgrade.

## Version 1.13.30 (2025-05-18)

- Update backend and CLI.

## Version 1.13.29 (2025-05-18)

- Update backend and CLI.

## Version 1.13.27 (2025-05-17)

- Update backend and CLI.

## Version 1.13.26 (2025-05-16)

- Update backend and CLI.

## Version 1.13.25 (2025-05-16)

- Update backend and CLI.

## Version 1.13.24 (2025-05-16)

- Update backend and CLI.

## Version 1.13.23 (2025-05-16)

- Update backend and CLI.

## Version 1.13.22 (2025-05-16)

- Update backend and CLI.

## Version 1.13.21 (2025-05-16)

- Update backend and CLI.

## Version 1.13.20 (2025-05-14)

- Update backend and CLI.

## Version 1.13.19 (2025-05-14)

- Update backend and CLI.

## Version 1.13.1 (2025-04-11)

- Update backend and CLI.

## Version 1.7.3 (2024-04-07)

- Preparing for next Bioconductor Release.

## Version 1.7.2 (2024-02-10)

- Fix vignette bug.

## Version 1.7.1 (2024-02-10)

- Update R version requirement.

## Version 1.5.2 (2023-08-24)

- Update java backend to BioInfoJavaUtils-1.4.0 (haploid GT in vcf).

## Version 1.5.1 (2023-04-30)

- Update java backend to BioInfoJavaUtils-1.3.1 (fasta header name).

## Version 1.1.6 (2022-09-26)

- Update java backend to BioInfoJavaUtils-1.2.4 (revamped FastaManager).

## Version 1.1.5 (2022-06-23)

- Update vignette minor bug.

## Version 1.1.4 (2022-06-16)

- Update vignette to handle getting sample files through https failure.

## Version 1.1.3 (2022-06-06)

- Update vignette to handle getting sample files through https failure.

## Version 1.1.2 (2022-05-15)

- Update possibly corrupted samples.vcf.gz.

## Version 1.1.1 (2022-05-08)

- Update tests to improve coverage.

## Version 1.0.0 (2022-04-27)

- Bioconductor 3.15 Release. New package **fastreeR**, Phylogenetic, Distance and Other Calculations on VCF and Fasta Files.

## Version 0.99.7 (2022-04-02)

- Updates and corrections after the 1st Bioconductor review.

## Version 0.99.6 (2022-03-26)

- Drop function `dist2hist`.

## Version 0.99.5 (2022-03-25)

- Update vignette's use of `dist2hist`.

## Version 0.99.4 (2022-03-25)

- Update java backend (createHistogram).

## Version 0.99.3 (2022-03-25)

- Update java backend.

## Version 0.99.2 (2022-03-25)

- Update `README.md` to inform about JDK>=8 requirement.
- Update java dependencies (jfreechart-1.5.3).

## Version 0.99.1 (2022-03-24)

- Update vignette to use `BiocFileCache` so that sample vcf and fasta downloads not get repeated needlessly.

## Version 0.99.0 (2022-03-21)

- Submitted to Bioconductor.
