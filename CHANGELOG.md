# Changelog

All notable changes to this project will be documented in this file.

This changelog is automatically updated by [release-please](https://github.com/googleapis/release-please) when contributors follow [conventional commit](https://www.conventionalcommits.org/) git messages. If you are using conventional commit messages, you should never need to edit this file manually. A github action will automatically update this file when a new release is created.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
