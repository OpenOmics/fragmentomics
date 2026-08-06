# Pipeline overview

## 1. About

This section documents what the pipeline *does* once <code>fragmentomics <b>run</b></code> has been invoked. Where the [Commands](../usage/run.md) section documents the command-line interface, this section documents the workflow behind it: how inputs are normalized, which reference files are used, and what each analysis step computes.

Every step runs inside a Singularity container defined in `config/containers.json`, and every step draws its threads, memory, walltime, partition and `--gres` request from `config/cluster.json`. No step depends on environment modules, so the workflow behaves identically under `--mode slurm` and `--mode local`.

## 2. Data flow

The pipeline accepts **either** paired-end Illumina FastQ files **or** ready-made alignments (BAM/CRAM/SAM), never both in the same run. The input type is auto-detected by the frontend and recorded in `config.json` as `project.input_type`. Both input types converge on a single canonical per-sample file:

```text
bams/{sample}.sorted.bam
```

Every analysis step reads that file and nothing else, which is what makes the analysis half of the pipeline input-type agnostic.

=== "FastQ input"

    ```text
    inputs/{sample}.R1.fastq.gz
    inputs/{sample}.R2.fastq.gz
             │
             │  align_fastq   (bwa-mem2 → fixmate → sort → markdup → filter)
             ▼
    bams/{sample}.sorted.bam
    ```

    See [FastQ alignment](fastq-alignment.md).

=== "BAM input"

    ```text
    <user-provided>.bam / .cram / .sam
             │
             │  stage_bams              (coordinate sort + index)
             ▼
    staged_bams/{sample}.sorted.bam
             │
             │  filter_reference_contigs  (subset to reference contigs)
             ▼
    bams/{sample}.sorted.bam
    ```

    See [BAM normalization](bam-normalization.md) and [Reference contig filter](contig-filter.md).

    !!! note

        `staged_bams/` only exists when the selected genome build ships a sequence dictionary (a `dict` entry in `config/genome.json`). Both bundled builds, `hg19` and `hg38`, do. If a build has no dictionary there is nothing to validate against, so `stage_bams` writes `bams/` directly and the contig filter is not part of the workflow at all.

From the canonical BAM, the analysis steps fan out in parallel:

```text
                       ┌─ coverage ───────────────► merge_coverage_excel
                       ├─ frag-length-bins
                       ├─ frag-length-intervals
                       ├─ end-motifs ─────────────► mds
bams/{sample}.sorted.bam ─┼─ interval-end-motifs
                       ├─ delfi
                       ├─ wps ──┬─────────────────► agg-bw  (aggregate)
                       │        └─ adjust-wps ────► agg-bw  (aggregate)
                       ├─ cleavage-profile ───────► agg-bw  (aggregate)
                       └─ bam_stats ──────────────► multiqc
```

Each analysis step is documented in the [Analyses](../analyses/index.md) section.

## 3. Which steps run

Some steps are conditional on the reference files the selected build provides. The workflow only requests their output if the files they need are present in `config/genome.json`, so a partially-populated build skips those analyses instead of failing on a missing path.

| Step | Condition |
|------|-----------|
| [`frag-length-intervals`](../analyses/frag-length-intervals.md) | a non-zero `--split-interval` |
| [`end-motifs`](../analyses/end-motifs.md), [`mds`](../analyses/mds.md) | `ref2bit` |
| [`interval-end-motifs`](../analyses/interval-end-motifs.md) | `ref2bit`, `intervals` |
| [`delfi`](../analyses/delfi.md) | `ref2bit`, `intervals`, `chrom_sizes` |
| [`wps`](../analyses/wps.md) | `tss` |
| [`agg-bw`](../analyses/agg-bw.md) on raw WPS | `tss`, `tss_interval` |
| [`adjust-wps`](../analyses/adjust-wps.md), [`cleavage-profile`](../analyses/cleavage-profile.md), and their aggregates | `tss`, `tss_interval`, `chrom_sizes` |
| [`filter_reference_contigs`](contig-filter.md) | BAM input **and** a `dict` entry |

The remaining steps are unconditional: [`frag-length-bins`](../analyses/frag-length-bins.md), [`coverage`](../analyses/coverage.md), and the project-level [QC steps](quality-control.md) (`bam_stats`, `merge_coverage_excel`, `multiqc`) always run. Note that `coverage` needs the build's `tss_interval` file even though it is not gated on it, so a build lacking that entry fails at the coverage step rather than skipping it.

Both bundled builds define every reference key, so a default `hg38` or `hg19` run executes every step.

## 4. Output directory layout

```text
<--output>/
├── inputs/                     symlinks to the user's input files (read-only)
├── staged_bams/                sorted BAMs awaiting the contig filter (BAM input only)
├── bams/                       canonical analysis BAM + index, one per sample
├── coverage/
│   ├── {sample}_coverage.bed
│   └── coverage_summary.xlsx   merged workbook across all samples
├── frag_length_bins/
│   ├── {sample}_frag_bin{bin_size}.tsv
│   └── {sample}_frag_bin{bin_size}.png
├── frag_length_intervals/{sample}_frag_interval.bed
├── end_motifs/{sample}_endmotif.tsv
├── interval_end_motifs/{sample}_endmotif_interval.tsv
├── mds/{sample}_mds.tsv
├── delfi/{sample}_delfi.bed
├── wps/
│   ├── {sample}_wps_out_tss.bw
│   └── {sample}_wps_out_tss_aggr.wig
├── adjust_wps/
│   ├── {sample}_wps_out_tss_adjusted.bw
│   └── {sample}_wps_out_tss_adj_aggr.wig
├── cleavage_profile/
│   ├── {sample}_cleavage_profile_tss.bw
│   └── {sample}_cleavage_profile_aggr.wig
├── qc/                         per-sample samtools reports + contig validation reports
├── multiqc/multiqc_report.html
├── config/                     resolved configuration for this run
├── workflow/                    Snakemake rules and scripts for this run
└── logfiles/                   master job log and per-job logs (slurm mode)
```

Sample names (`{sample}`) are derived from the input file basenames with the `.R1/.R2` mate suffix or the `.sorted`/`.bam`/`.cram`/`.sam` extensions stripped. Two inputs that reduce to the same sample name are rejected up front rather than silently clobbering one another.

## 5. Reproducibility

The `--output` directory is self-contained. `config/`, `workflow/` and the resolved `config.json` are copied into it at launch, so a completed output directory records exactly which rules, scripts, reference paths and resource requests produced its results. Re-running in the same directory resumes rather than restarts; passing `--overwrite-pipeline-template` refreshes the copied template from the current installation.

See also: [Reference files](references.md) and [Quality control](quality-control.md).
