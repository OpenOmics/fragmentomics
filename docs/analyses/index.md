# Analyses

## 1. About

The analysis half of the pipeline runs [FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit), a Python package for cfDNA fragmentation analysis. Each page in this section documents one FinaleToolkit subcommand as this pipeline invokes it: what it measures, exactly which flags the pipeline passes, what the output looks like, and which `--genome` reference files and <code>fragmentomics <b>run</b></code> options control it.

Every analysis reads the same input — the canonical per-sample analysis BAM, `bams/{sample}.sorted.bam` — regardless of whether that BAM came from [aligning FastQ files](../pipeline/fastq-alignment.md) or from [staging and filtering user BAMs](../pipeline/bam-normalization.md).

!!! note "These are workflow steps, not commands you run"

    You do not invoke these subcommands yourself. <code>fragmentomics <b>run</b></code> constructs and submits them. The flag listings on these pages document what the pipeline runs on your behalf, so you can interpret the outputs and know which `run` options change them.

## 2. The analyses

<section markdown="1">

| Analysis | Measures | Output |
|----------|----------|--------|
| [`coverage`](coverage.md) | Normalized fragment coverage per interval | `coverage/{sample}_coverage.bed` |
| [`frag-length-bins`](frag-length-bins.md) | Genome-wide fragment-length histogram | `frag_length_bins/{sample}_frag_bin{n}.tsv` + `.png` |
| [`frag-length-intervals`](frag-length-intervals.md) | Fragment-length statistics per interval | `frag_length_intervals/{sample}_frag_interval.bed` |
| [`end-motifs`](end-motifs.md) | Genome-wide 4-mer fragment-end motif frequencies | `end_motifs/{sample}_endmotif.tsv` |
| [`interval-end-motifs`](interval-end-motifs.md) | End-motif frequencies per interval | `interval_end_motifs/{sample}_endmotif_interval.tsv` |
| [`mds`](mds.md) | Motif diversity score (one number per sample) | `mds/{sample}_mds.tsv` |
| [`delfi`](delfi.md) | GC-corrected short-to-long fragment ratios | `delfi/{sample}_delfi.bed` |
| [`wps`](wps.md) | Windowed protection score around TSSs | `wps/{sample}_wps_out_tss.bw` |
| [`adjust-wps`](adjust-wps.md) | Smoothed, edge-corrected WPS | `adjust_wps/{sample}_wps_out_tss_adjusted.bw` |
| [`cleavage-profile`](cleavage-profile.md) | Per-base cleavage proportion around TSSs | `cleavage_profile/{sample}_cleavage_profile_tss.bw` |
| [`agg-bw`](agg-bw.md) | Aggregate of a bigWig signal across all windows | `*_aggr.wig` |

</section>

## 3. Shared conventions

Most of the analyses take the same handful of options, and the pipeline sets them consistently:

  `-q 30` — mapping quality threshold
> **Applied by every analysis that reads the BAM.**
>
> Fragments whose reads map below MAPQ 30 are excluded. This is applied at read time rather than baked into the BAM, so the same analysis BAM can be re-analyzed at a different threshold without realigning.

---
  `-min` / `-max` — fragment length window
> **Set from `--fragment-minimum` and `--fragment-maximum`** (defaults 50 and 500).
>
> Applied by `coverage`, `frag-length-bins`, `frag-length-intervals`, `end-motifs`, `interval-end-motifs` and `cleavage-profile`. [`wps`](wps.md) is the exception: it uses its own fixed 120–180 bp window, which is what makes it an *L*-WPS.

---
  `-p` — intersect policy
> **How a fragment is assigned to an interval.**
>
> `any` (overlap anywhere) for [`coverage`](coverage.md) and [`frag-length-intervals`](frag-length-intervals.md); `midpoint` (the fragment's midpoint must fall inside) for [`frag-length-bins`](frag-length-bins.md).

---
  `-w` — worker processes
> **Set from the step's `threads` in `config/cluster.json`.**
>
> Every analysis that supports parallel workers gets its thread allocation passed through, so tuning a step's throughput is a matter of editing `config/cluster.json` in the output directory.

---
  `-v` — verbose
> **Always enabled.**
>
> Progress and parameter detail go to the step's log, which is where to look when an analysis produces unexpectedly empty output.

## 4. Other FinaleToolkit subcommands

FinaleToolkit 0.11.0 also provides `delfi-gc-correct`, `breakpoint-motifs`, `interval-breakpoint-motifs`, `interval-mds`, `filter-file` and `gap-bed`. This pipeline does not currently run them, so they are not documented here — see the [FinaleToolkit documentation](https://finaletoolkit.readthedocs.io/) if you need them.
