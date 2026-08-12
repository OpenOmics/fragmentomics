# Quality control

## 1. About

The pipeline produces several kinds of quality-control output. The `bam_stats` reports and the project-level outputs are written regardless of which input path was taken; the FastQC, `fastp` and duplicate reports exist only for a FastQ run, which is the only path that has raw reads to inspect and trim and that does its own duplicate marking:

| Output | Produced by | Scope |
|--------|-------------|-------|
| `qc/{sample}.samtools.stats.txt`<br>`qc/{sample}.flagstat.txt`<br>`qc/{sample}.idxstats.txt` | `bam_stats` | per sample |
| `qc/{sample}.R1_fastqc.{zip,html}`<br>`qc/{sample}.R2_fastqc.{zip,html}`<br>`qc/{sample}.sorted_fastqc.{zip,html}`<br>`qc/{sample}.markdup.stats.txt`<br>`qc/{sample}.fastp.{json,html}` | [`align_fastq`](fastq-alignment.md) | per sample, FastQ input only |
| `qc/{sample}.contig_validation.txt` | [`filter_reference_contigs`](contig-filter.md) | per sample, BAM input only |
| `multiqc/multiqc_report.html` | `multiqc` | whole project |
| `coverage/coverage_summary.xlsx` | `merge_coverage_excel` | whole project |

## 2. Alignment metrics (`bam_stats`)

`bam_stats` runs `samtools` against the canonical analysis BAM for each sample and writes three reports:

  `samtools stats`
> **Comprehensive alignment summary.**
>
> Read counts, mapping rates, error rate, insert-size distribution, GC content, coverage distribution, and per-cycle base composition. This is the richest of the three and supplies most of the MultiQC plots.

---
  `samtools flagstat`
> **Flag-level tallies.**
>
> Counts by SAM flag: total, mapped, properly paired, duplicates, singletons, and mate-on-different-contig. This is where the **duplicate rate** comes from — see §2.1.

---
  `samtools idxstats`
> **Per-contig read counts.**
>
> One line per contig with its length and its mapped and unmapped read counts. Useful for spotting a contig that is unexpectedly empty, and for confirming after a BAM-input run that only the reference's contigs remain.

### 2.1 Why samtools and not Picard

Metrics are derived from `samtools` rather than from Picard because the [FastQ path](fastq-alignment.md) fuses alignment, filtering, fixmate, sorting and duplicate marking into a **single pipe**, and `samtools markdup` reports on that marking directly with `-f` (see §2.2). `flagstat` additionally recovers the duplicate rate from the flags on the reads themselves, which works identically for both input paths — including for staged BAMs, whose duplicate flags were set by whatever tool the user ran before providing them.

All three reports are natively parsed by MultiQC, so the aggregate report has content no matter which path produced the BAM.

## 2.2 FastQ-path reports (`align_fastq`)

A FastQ run adds three more kinds of report, all written into `qc/` by [`align_fastq`](fastq-alignment.md) and all parsed natively by MultiQC:

  `samtools markdup -f`
> **Duplicate report,** `qc/{sample}.markdup.stats.txt`.
>
> Reads read/written/excluded/examined, paired and single duplicate counts, optical duplicates, and an estimated library size. Because `markdup` is the last stage of the pipe, these duplicates are flagged in the analysis BAM rather than removed from it, so this report describes reads that downstream analyses will still see.

---
  `fastqc`
> **Read QC before and after alignment.**
>
> `qc/{sample}.R1_fastqc.*` and `qc/{sample}.R2_fastqc.*` from the raw mates, and `qc/{sample}.sorted_fastqc.*` from the analysis BAM. The first pair characterizes the library as sequenced — untrimmed, since FastQC runs upstream of `fastp`; the third characterizes what survived trimming, filtering and alignment. Adapter content is the clearest before/after contrast.

---
  `fastp`
> **Adapter-trimming and base-quality filtering report,** `qc/{sample}.fastp.json` (plus an HTML copy).
>
> How many reads had adapter read-through trimmed and how many bases that removed, plus how many read pairs `fastp` kept and how many it dropped for falling below the [`--baseqscore`](../usage/run.md#22-analysis-options) mean-quality threshold or for trimming shorter than 15 bp — with before/after read counts, base counts, read lengths, Q20/Q30 rates and GC content. MultiQC identifies the JSON by content, not by filename. This is the report that accounts for the gap between the input FastQ read count and the reads `bwa-mem2` ever saw — the [`samtools view`](fastq-alignment.md#34-adapter-trimming-and-read-filtering) filter and duplicate marking account for the rest.

## 3. Aggregate report (`multiqc`)

`multiqc` scans the `qc/` directory and writes `multiqc/multiqc_report.html` plus its `multiqc_report_data/` directory. Only `qc/` is scanned — deliberately, so MultiQC does not walk the large bigWig and BED outputs of the analysis steps looking for something to parse.

The report is a project-level gather: it waits for every sample's `bam_stats` reports, for the FastQC, `fastp` and markdup reports on a FastQ run, and for the merged coverage workbook, so it is written once per pipeline invocation with all samples side by side. The contig-validation reports share the `qc/` directory but are not in a format MultiQC recognizes, so they are ignored.

## 4. Merged coverage workbook (`merge_coverage_excel`)

Every per-sample [`coverage`](../analyses/coverage.md) BED is merged into a single Excel workbook, `coverage/coverage_summary.xlsx`, with two sheets:

  `summary`
> **One row per sample.**
>
> Interval count, total coverage, and the mean, median, standard deviation, minimum and maximum of the per-interval coverage values. This is the sheet to check first for a sample that is an outlier relative to the rest of the cohort.

---
  `coverage`
> **The interval × sample matrix.**
>
> One row per genomic interval, one column per sample, holding the normalized coverage value. Intervals are the outer join across all samples, so an interval present in only one sample appears with blanks elsewhere. Rows are ordered by first-appearance contig order and then by coordinate, so `chr1` precedes `chr2` precedes `chr10` rather than sorting lexically.

Two behaviors worth knowing:

- **Duplicate sample names are rejected.** Two BEDs reducing to the same sample name would silently overwrite one another's column, so the step fails instead.
- **Large interval sets fall back to TSV.** An Excel worksheet holds at most 1,048,576 rows. If the merged matrix exceeds that, the workbook is still written but the `coverage` sheet is emitted as a `.tsv` sidecar alongside it rather than being truncated.

## 5. What to check after a run

  **Did the right reference get used?** *(BAM input)*
> Open `qc/{sample}.contig_validation.txt`. `length_mismatches` must be `0` — the run would have failed otherwise — and `matched` should equal `reference_contigs` for a normal whole-genome BAM. See [reading the outcome](contig-filter.md#7-reading-the-outcome).

---
  **Is the library complex enough?**
> Duplicate rate in the MultiQC general statistics table, from `flagstat` and — on a FastQ run — from `qc/{sample}.markdup.stats.txt`. A high rate means fragmentation features are being computed from fewer independent molecules than the raw read count suggests. On a FastQ run the duplicates are marked but retained in the analysis BAM, so they contribute to coverage and to the fragment-length and motif distributions.

---
  **Are fragment lengths as expected?**
> The insert-size plot in MultiQC, and per sample the histogram written by [`frag-length-bins`](../analyses/frag-length-bins.md) at `frag_length_bins/{sample}_frag_bin{bin_size}.png`. cfDNA should show the characteristic ~167 bp mononucleosome peak with a shoulder near 320 bp. A sample lacking it is likely contaminated with genomic DNA.

---
  **Is coverage comparable across samples?**
> The `summary` sheet of `coverage/coverage_summary.xlsx`. Coverage is normalized and scaled, so mean values should be broadly comparable; a sample far off the others warrants investigation before its features are compared to the rest.

---
  **Did any contig come out empty?**
> `qc/{sample}.idxstats.txt`. A primary chromosome with zero mapped reads points at a reference or input problem rather than biology.

---
  **How much of the data did filtering remove?**
> On a FastQ run, `qc/{sample}.fastp.json` gives the pairs dropped before alignment — for mean base quality, or for trimming below 15 bp — and the gap between that and `flagstat`'s total gives what the `samtools view` flag and MAPQ filter removed after it. Losing a large fraction to either is a signal to revisit [`--baseqscore` and `--mapscore`](../usage/run.md#22-analysis-options) rather than to accept a thin BAM; on a BAM run the same thresholds are applied during [staging](bam-normalization.md#22-read-filtering), where the only record of them is the step's log.
