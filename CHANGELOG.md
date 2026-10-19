# Changelog

All notable changes to this project will be documented in this file.

This changelog is automatically updated by [release-please](https://github.com/googleapis/release-please) when contributors follow [conventional commit](https://www.conventionalcommits.org/) git messages. If you are using conventional commit messages, you should never need to edit this file manually. A github action will automatically update this file when a new release is created.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0](https://github.com/OpenOmics/fragmentomics/compare/v1.0.0...v1.1.0) (2026-10-19)


### Features

* **CDA:** dockerfile for cfDNAanalyzer ([25fcfef](https://github.com/OpenOmics/fragmentomics/commit/25fcfef4e6f42408491ca33fb8ee747509c1a721))
* **CDA:** snakemake routing ([1effa87](https://github.com/OpenOmics/fragmentomics/commit/1effa8739f84070de9c1bdff1656694fec4bb849))
* **cfdna:** documentation ([7e1ad61](https://github.com/OpenOmics/fragmentomics/commit/7e1ad61df680f31d9f05475c461e43df602854fe))
* **cfdna:** implement cfdna feature extraction pipeline ([7261a94](https://github.com/OpenOmics/fragmentomics/commit/7261a94373873694174e9f4ebdb49974473e4ad7))
* **cfdna:** infrastructure scripting ([90d7882](https://github.com/OpenOmics/fragmentomics/commit/90d7882d6b8d14e9d5e38a289585b9e9e09960ad))
* **docker:** add dockerfile for cfdnaanalyzer ([6d59168](https://github.com/OpenOmics/fragmentomics/commit/6d5916890412a2bcc20c30e544b455266cc37946))


### Bug Fixes

* **cfdnaanalyzer:** ignore header-only inputs during schema construction ([a9c3b4d](https://github.com/OpenOmics/fragmentomics/commit/a9c3b4dfb64b57a7f6a0a1dff7acbf8916bdbf3b))
* **cfdnaanalyzer:** use bounded-memory matrix joins ([c75cfb5](https://github.com/OpenOmics/fragmentomics/commit/c75cfb5d8330de76361129612534308c517cd455))
* **docker:** detach R installs ([240afef](https://github.com/OpenOmics/fragmentomics/commit/240afef747fbaaf04571c559752ae9f159608eb0))
* **docker:** inline R library installs ([cb0b59a](https://github.com/OpenOmics/fragmentomics/commit/cb0b59a7b94469c0d65d6ac7d3fbb130e4d241bf))
* **docker:** ncurses bc and libtinfo5 ([89904fa](https://github.com/OpenOmics/fragmentomics/commit/89904faa84977ea05c9ddca839500de440b780ea))
* **docker:** r2py configuration ([0780154](https://github.com/OpenOmics/fragmentomics/commit/0780154d8f98bba631be56e97b06c1d1d735cb92))
* **run.sh:** revert job time change ([adbc735](https://github.com/OpenOmics/fragmentomics/commit/adbc735e295dfd62dffda9a0bc6db0a78ca837bb))

## [1.0.0](https://github.com/OpenOmics/fragmentomics/releases/tag/v1.0.0) (2026-08-25)

### ⚠ BREAKING CHANGES

* **input:** SAM and CRAM inputs are no longer supported. The pipeline now accepts paired-end Illumina FastQ files or BAM alignments ([#22](https://github.com/OpenOmics/fragmentomics/pull/22)).

### Features

* **alignment:** add paired-end FastQ processing with adapter trimming, quality filtering, `bwa-mem2` alignment, mate repair, coordinate sorting, duplicate marking, indexing, and FastQC reporting ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **alignment:** add CPU-aware `bwa-mem2` binary selection for Intel and AMD compute nodes ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **input:** add unified FastQ and BAM entry points that produce a canonical analysis-ready BAM ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **references:** add complete `hg19` and `hg38` reference configurations, including alignment indexes, FASTA files, sequence dictionaries, chromosome sizes, TSS annotations, blacklist regions, and gap tracks ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **bam:** add BAM normalization, orphaned-mate removal, and reference-contig validation and filtering ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **cli:** add `--mapscore`, `--baseqscore`, and `--interval` options ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **intervals:** generate fixed-width genomic intervals at runtime with support for SI-prefixed sizes such as `5kb`, `1mb`, and `2MB` ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **output:** add BGZF-compressed and tabix-indexed BED exports of analysis alignments ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **qc:** add samtools, FastQC, fastp, and MultiQC quality-control reporting ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **multiqc:** add custom fragmentomics sections for fragment lengths, end motifs, MDS, DELFI, coverage, WPS, and cleavage profiles ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **coverage:** add a merged Excel workbook containing per-sample statistics and an interval-by-sample coverage matrix ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **docker:** add dedicated `bwa-mem2` and FinaleToolkit containers with pinned workflow dependencies ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))

### Bug Fixes

* **intervals:** prevent runs with different window sizes from reusing or overwriting incompatible interval files ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **local:** resolve SLURM-based temporary-directory settings when running locally ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* **validation:** improve sample-name collision detection during staging and coverage aggregation ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))

### Documentation

* document pipeline architecture, alignment, BAM normalization, references, quality control, output layout, and individual fragmentomics analyses ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
* update command-line documentation for FinaleToolkit 1.1.0 ([#22](https://github.com/OpenOmics/fragmentomics/pull/22))
