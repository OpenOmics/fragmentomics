# Quality control

## 1. About

The pipeline produces several kinds of quality-control output. Both input paths produce alignment metrics, a FastQC report on the analysis BAM and a `fastp` report, but from different rules: on a FastQ run they fall out of the alignment work itself, while on a BAM run they are produced by two rules that exist for no other purpose. Only a FastQ run has raw reads to inspect, and only a FastQ run marks duplicates itself, so the pre-alignment FastQC reports and the `markdup` report are specific to it:

| Output | Produced by | Scope |
|--------|-------------|-------|
| `qc/{sample}.samtools.stats.txt`<br>`qc/{sample}.flagstat.txt`<br>`qc/{sample}.idxstats.txt` | `bam_stats` | per sample |
| `qc/{sample}.R1_fastqc.{zip,html}`<br>`qc/{sample}.R2_fastqc.{zip,html}`<br>`qc/{sample}.sorted_fastqc.{zip,html}`<br>`qc/{sample}.markdup.stats.txt`<br>`qc/{sample}.fastp.{json,html}` | [`align_fastq`](fastq-alignment.md) | per sample, FastQ input only |
| `qc/{sample}.sorted_fastqc.{zip,html}` | `fastqc_bam` | per sample, BAM input only |
| `qc/{sample}.fastp.{json,html}` | `fastp_bam` | per sample, BAM input only |
| `qc/{sample}.contig_validation.txt` | [`filter_reference_contigs`](contig-filter.md) | per sample, BAM input only |
| `qc/finaletoolkit/` | `finaletoolkit_multiqc` | whole project |
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

## 2.3 BAM-path reports (`fastqc_bam`, `fastp_bam`)

A BAM run has no raw reads, and its filtering happened during [staging](bam-normalization.md#22-read-filtering) with `samtools` rather than with `fastp`, so neither of the reports above comes for free. Both are produced anyway, from the analysis BAM, by two rules whose only purpose is reporting — nothing they do changes the BAM. They write to the same filenames as their FastQ-path counterparts, so MultiQC labels the sample identically and the same sections mean the same thing on either path:

  `fastqc --format bam`
> **Read QC of the analysis BAM,** `qc/{sample}.sorted_fastqc.*`.
>
> The same report `align_fastq` produces after alignment, on the same file, since FastQC reads BAM natively: per-base and per-sequence quality, GC content, read-length distribution, adapter and overrepresented-sequence content. What a FastQ run has in addition is the *before* half of that comparison, from the raw mates.
>
> One caveat applies to this report on both paths: FastQC estimates duplication and overrepresented sequences from the first 100,000 reads, which in a coordinate-sorted BAM all come from the start of the first contig rather than from across the library. Read the quality and content modules here, and take the duplicate rate from `flagstat`.

---
  `fastp` on the BAM converted back to FastQ
> **Descriptive read report,** `qc/{sample}.fastp.json` (plus an HTML copy).
>
> The BAM is streamed back to FastQ (`samtools collate` → `samtools fastq`, name-collated so the mate files line up, staged in the node's temporary directory and discarded afterwards) and `fastp` is pointed at the result with **every transformation disabled** — no adapter trimming, no quality filter, no length filter, no polyG trimming. It only measures, and writes no reads back out.
>
> Everything descriptive is therefore meaningful: per-base quality and content curves, Q20/Q30 rates, GC content, duplication rate, and the insert-size distribution `fastp` estimates from mate overlap. What is *not* meaningful is the before/after comparison — the two halves are identical by construction, and the trimming and filtering counts are all zero. That is deliberate: with trimming left on, the report would show adapter bases removed and reads dropped that no step of a BAM run actually loses, and the same MultiQC section would then mean the opposite of what it means on a FastQ run, where those removals are real. Adapter content for a BAM run is in the FastQC report above instead. For the same reason the general statistics table has no `% Adapter` value for a BAM run.
>
> Reads whose mate did not survive filtering are left out, since they cannot be reported as pairs; there are normally very few. A BAM with no paired reads at all is reported single-end rather than failing the step — the layout is read from the sample's `flagstat` report.

## 3. Aggregate report (`multiqc`)

`multiqc` scans the `qc/` directory and writes `multiqc/multiqc_report.html` plus its `multiqc_report_data/` directory. Only `qc/` is scanned — deliberately, so MultiQC does not walk the large bigWig and BED outputs of the analysis steps looking for something to parse.

The report is a project-level gather: it waits for every sample's `bam_stats` reports, for the read QC reports of whichever input path ran (§2.2, §2.3), for the fragmentomics summary (§3.1), and for the merged coverage workbook, so it is written once per pipeline invocation with all samples side by side. The contig-validation reports share the `qc/` directory but are not in a format MultiQC recognizes, so they are ignored.

### 3.1 Fragmentomics sections (`finaletoolkit_multiqc`)

MultiQC has no `finaletoolkit` module, so nothing the analysis steps produce would reach the report on its own — their tables, wigs and bigWigs are simply files MultiQC does not recognize. `finaletoolkit_multiqc` reads the per-sample results and rewrites what is summarizable as MultiQC *custom content* into `qc/finaletoolkit/`, one document per section, which the scan of `qc/` then picks up. Nothing is recomputed and no BAM is re-read — it is a summary of files that already exist, and the largest of them, the per-interval end motif tables, are read in chunks rather than held whole.

Four headline numbers join the **general statistics** table next to the samtools and FastQC columns — median fragment length, the percentage of fragments under 150 bp, [MDS](../analyses/mds.md) and mean [coverage](../analyses/coverage.md) — and the rest is grouped under a **Fragmentomics** heading:

| Section | Derived from |
|---------|--------------|
| Fragmentation metrics (table) | all of the below, one row per sample |
| Fragment length distribution | [`frag-length-bins`](../analyses/frag-length-bins.md) |
| Fragment length across intervals | [`frag-length-intervals`](../analyses/frag-length-intervals.md) |
| End motif frequency (top 12) | [`end-motifs`](../analyses/end-motifs.md) |
| End motif frequency across intervals (top 12) | [`interval-end-motifs`](../analyses/interval-end-motifs.md) |
| DELFI fragmentation profile | [`delfi`](../analyses/delfi.md) |
| TSS WPS profile | [`wps`](../analyses/wps.md) aggregate |
| TSS adjusted WPS profile | [`adjust-wps`](../analyses/adjust-wps.md) aggregate |
| TSS cleavage profile | [`cleavage-profile`](../analyses/cleavage-profile.md) aggregate |

Every section is an overlay of all samples in the run, and each carries its own note on how to read it. Two of them summarize *intervals* rather than fragments and are worth a word here: **Fragment length across intervals** is the distribution of the per-interval median length, so it says how uniform the fragmentation is along the genome rather than how long the fragments are, and **End motif frequency across intervals** is the motif profile of the fragments that fell in the target intervals, weighted by the fragment ends in each interval — read against the genome-wide section above it, a difference between the two is a difference between the intervals and the rest of the genome.

What is **not** in the report is anything per-interval, per-bin or per-base: the interval × sample coverage matrix (that is the [workbook](#4-merged-coverage-workbook-merge_coverage_excel)), the per-interval tables, the BED of aligned intervals, and the per-sample bigWigs, which are for a genome browser rather than for a report. The sections above are summaries of those files, not replacements for them.

Which sections exist depends on the run: most analysis steps only run when the selected genome build supplies the reference files they need (see [Reference files](references.md)), and the interval fragment lengths additionally need a non-zero [`--split-interval`](../usage/run.md#22-analysis-options). The summary covers whichever steps did run. Fragment length and coverage need nothing beyond the analysis BAM, so those are always present.

The profile plots are downsampled to keep the report interactive: the DELFI profile has hundreds of thousands of 5 kb bins per sample, averaged into 500 points of a few megabases each — the scale the fragmentation profile is normally read at — and the TSS profiles are reduced to roughly one point per 8 bp of their ±2 kb window. The full-resolution values are in the analysis outputs themselves; these sections are for comparing samples at a glance, not for measuring.

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
>
> Take the rate from those two and not from the FastQC or `fastp` columns beside them: FastQC estimates it from the first 100,000 reads of the file, and `fastp` from identical sequences rather than from alignment positions, so both differ from the flag-based rate — by design, not because something is wrong.

---
  **Are fragment lengths as expected?**
> The **Fragment length distribution** section of the MultiQC report (§3.1), which overlays every sample, with the insert-size plots from `samtools` and `fastp` as a cross-check; per sample, the histogram written by [`frag-length-bins`](../analyses/frag-length-bins.md) at `frag_length_bins/{sample}_frag_bin{bin_size}.png`. cfDNA should show the characteristic ~167 bp mononucleosome peak with a shoulder near 320 bp. A sample lacking it is likely contaminated with genomic DNA. The median fragment length and the percentage of fragments under 150 bp are in the general statistics table for a quick scan across the cohort.

---
  **Is coverage comparable across samples?**
> The `summary` sheet of `coverage/coverage_summary.xlsx`. Coverage is normalized and scaled, so mean values should be broadly comparable; a sample far off the others warrants investigation before its features are compared to the rest.

---
  **Did every read keep its mate?**
> *Singletons* in `qc/{sample}.flagstat.txt` should be `0` on both input paths — the FastQ path guarantees it with the `fixmate` guard in [§3.5 of FastQ alignment](fastq-alignment.md#35-orphaned-mates), and the BAM path with [`remove_orphan_reads`](bam-normalization.md#3-what-remove_orphan_reads-does). Every fragmentation feature here is computed from a pair's coordinates, so a surviving orphan is a fragment inferred from a mate that is not in the file.
>
> A non-zero count on a BAM run means the input was single-end and was passed through unchanged rather than emptied; the step logs a warning saying so, and the fragment-level analyses have nothing meaningful to measure on such a BAM. Note that unlike the FastQ path, *properly paired* is **not** expected to equal *total* on a BAM run — discordant pairs are deliberately kept (see [what is not done](bam-normalization.md#4-what-is-not-done)).

---
  **Did any contig come out empty?**
> `qc/{sample}.idxstats.txt`. A primary chromosome with zero mapped reads points at a reference or input problem rather than biology.

---
  **Do the fragmentomics features agree across samples?**
> The **Fragmentomics** sections of the MultiQC report (§3.1). The end-motif and TSS profiles should have the same shape for every sample of the same assay — a sample whose profile is flat where the others are structured, or whose MDS sits apart from the cohort, is worth investigating before its features are compared to the rest. Read the profiles as overlays for shape, not for exact values: two of them are downsampled for the report.

---
  **How much of the data did filtering remove?**
> On a FastQ run, `qc/{sample}.fastp.json` gives the pairs dropped before alignment — for mean base quality, or for trimming below 15 bp — and the gap between that and `flagstat`'s total gives what the `samtools view` flag and MAPQ filter removed after it. Losing a large fraction to either is a signal to revisit [`--baseqscore` and `--mapscore`](../usage/run.md#22-analysis-options) rather than to accept a thin BAM; on a BAM run the same thresholds are applied during [staging](bam-normalization.md#22-read-filtering), where the only record of them is the step's log. A BAM run then loses a further, usually small, share to the [pairing cleanup](bam-normalization.md#3-what-remove_orphan_reads-does) — the mates the per-record filters orphaned, plus whatever singletons the input arrived with — and the contig filter reports its own losses separately in `qc/{sample}.contig_validation.txt`.
