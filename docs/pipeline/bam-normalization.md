# BAM normalization

## 1. About

When the pipeline is given ready-made alignments instead of FastQ files, it runs the *full pipeline minus alignment*. The provided files still need to be brought into a uniform shape before the analysis steps can read them, because FinaleToolkit requires a coordinate-sorted, indexed BAM and user-supplied alignments arrive in whatever form their producer left them in.

Normalization is a two-stage process:

```text
<user-provided>.bam / .cram / .sam
         │
         │  stage_bams                 ← this page
         ▼
staged_bams/{sample}.sorted.bam
         │
         │  filter_reference_contigs   ← see "Reference contig filter"
         ▼
bams/{sample}.sorted.bam
```

This page covers the first stage. The second is documented in [Reference contig filter](contig-filter.md).

## 2. What `stage_bams` does

`stage_bams` runs once per pipeline invocation and processes every input file. For each one it:

1. **Determines the format from the file extension** — `.bam`, `.cram` or `.sam`, case-insensitively. A file whose extension is none of these is a hard error; the pipeline does not attempt to sniff the format from the file contents.
2. **Normalizes the sample name.** The basename is stripped of its `.sorted` marker (if present) and of its `.bam`/`.cram`/`.sam` extension, then `.sorted.bam` is appended. So `Sample_A.cram`, `Sample_A.sam` and `Sample_A.sorted.bam` all normalize to the sample name `Sample_A`.
3. **Filters reads on mapping and base quality** with `samtools view`, using the thresholds given on the command line. See [§2.2](#22-read-filtering).
4. **Coordinate sorts and converts to BAM** with `samtools sort`. CRAM and SAM inputs are converted to BAM as part of this step — there is no separate conversion step. Already coordinate-sorted input still passes through the sort, which is cheap relative to re-sorting and guarantees the ordering rather than trusting the header's `SO:` tag. The filter above is piped straight into the sort, so no intermediate alignment file is written to disk.
5. **Indexes the result**, producing a `.bai` alongside each BAM.
6. **Refreshes the index timestamp** so the index is never older than the BAM it indexes. Without this, `samtools` and `pysam` emit "index is older than data file" warnings on some filesystems where the two writes land in the same timestamp granularity.

The output lands in `staged_bams/`.

### 2.1 Duplicate sample names are rejected

Because sample names are derived from basenames, two inputs from different directories can collide:

```text
/path/one/Sample_A.bam   ─┐
                          ├─► both normalize to Sample_A
/path/two/Sample_A.cram  ─┘
```

`stage_bams` detects this before writing anything and fails with the list of colliding names. It does not silently pick one, and it does not disambiguate by adding a suffix — a collision means the two files would have produced results attributed to the same sample, and only the user knows which name each should have. Rename the inputs and re-run.

### 2.2 Read filtering

Two quality thresholds are applied while staging, both taken from the command line — [`--mapscore`](../usage/run.md#22-analysis-options) and [`--baseqscore`](../usage/run.md#22-analysis-options). Both are inclusive lower bounds, so a read is kept when its score is greater than or equal to the value given, and passing `0` for either leaves that filter out of the `samtools view` invocation entirely:

```text
samtools view -b -q <mapscore> -e 'avg(qual) >= <baseqscore>' <input> \
  | samtools sort -o staged_bams/<sample>.sorted.bam -
```

  `-q <mapscore>`
> **Drop alignments below the mapping-quality threshold.** Default 20. This is the same filter the [FastQ path](fastq-alignment.md#34-adapter-trimming-and-read-filtering) applies after alignment, so both input paths reach the analysis steps with the same MAPQ floor. The same value is also passed to each FinaleToolkit analysis as its `-q` threshold, so one setting governs the whole run.

---
  `-e 'avg(qual) >= <baseqscore>'`
> **Drop reads whose mean base quality is below the threshold.** Default 20. The FastQ path enforces the identical criterion with `fastp --average_qual` before alignment; here it has to be a filter expression, because by the time the pipeline sees a BAM the reads are already aligned. Note that a record carrying no quality string (`*`) is **kept** regardless of the threshold: htslib stores a missing quality as `0xff` per base, so `avg(qual)` evaluates to 255 for such a record and it passes any threshold. A BAM stripped of its quality strings is therefore not filtered on base quality at all.

Both filters are applied by `samtools` rather than `pysam.sort`, so the staging step requires `samtools` on `PATH`; it fails immediately with a clear message if filtering was requested and the binary is missing. When both thresholds are `0` the step falls back to `pysam.sort` and no filtering process is started at all.

!!! warning "Per-read filtering on the BAM path can orphan mates"

    These filters act on individual records, so a pair whose R1 passes and whose R2 fails leaves a record in the staged BAM whose mate is gone, and whose `0x2` (properly-paired) flag still claims otherwise. FinaleToolkit infers a fragment from the pair's coordinates, so such records contribute no fragment and are effectively discarded downstream — but they do still count in read totals in the QC reports.

    The FastQ path does not have this problem. `fastp` filters *pairs*, so both mates leave together, and the alignment step additionally runs [`samtools fixmate`](fastq-alignment.md#35-orphaned-mates) after its own filter to strip and drop any mate left behind. That guard cannot be applied here: `fixmate` needs name-collated input, and re-collating a staged BAM would mean a second whole-file sort on every input. If mate integrity in the staged BAM matters to you, filter and re-pair the input yourself and run with `--mapscore 0 --baseqscore 0`.

## 3. What is *not* done

Apart from the quality filtering in [§2.2](#22-read-filtering), `stage_bams` deliberately does not modify the read content of the input:

- **Duplicates are not marked.** The [FastQ path](fastq-alignment.md) runs `samtools markdup` so duplicates are at least *flagged* (it does not remove them either); staged BAMs are taken as-is, with whatever `0x400` flags — or none — they arrived with. If your BAMs have not been deduplicated, coverage will be inflated and fragment-length and end-motif distributions will be biased toward whatever the PCR amplified. Deduplicate before providing them.
- **No flag-based filtering is applied.** The FastQ path additionally keeps only properly-paired, primary, mapped reads (`-f 2 -F 2828`) and drops mates orphaned by its own filter. Staged BAMs keep every record that passes the quality thresholds, whatever its flags. FinaleToolkit re-applies the mapping-quality threshold and its fragment-length window at read time, but it does not remove supplementary or secondary alignments for you.
- **Read groups are not added or rewritten.** Whatever `@RG` lines the input carries are preserved. So are `@PG` lines, which means the provenance of the original alignment survives into the analysis BAM.
- **Reads are not re-aligned.** The alignment coordinates in a staged BAM are exactly the ones its original aligner produced, against whatever reference that aligner used. The pipeline cannot detect a *wrong* reference from the coordinates alone — that is what the [contig filter](contig-filter.md) is for.

!!! warning "Staged BAMs are assumed to be aligned to the selected `--genome` build"

    `--genome` is required for BAM input, because it selects the interval, TSS, chrom.sizes, 2bit and blacklist files used by the analysis steps. If your BAMs were aligned to a different build, those reference files describe different coordinates than your reads do and the results will be meaningless. The [contig filter](contig-filter.md) catches the common, detectable version of this mistake — a contig present in both with a different length — and fails the run. It cannot catch a build difference that leaves all shared contig lengths identical.

## 4. Resources

```json
"stage_bams": {
    "threads": 48,
    "mem": "96G",
    "time": "12:00:00",
    "partition": "norm"
}
```

`stage_bams` is a single step that loops over all inputs serially, so its walltime scales with the *total* size of your input set rather than the largest single file. For large cohorts, raise `time` in `config/cluster.json` in the output directory. The threads are passed to both `samtools view` and `samtools sort`, which run concurrently as two ends of a pipe.

## 5. Next step

Once staging completes, `staged_bams/{sample}.sorted.bam` is compared against the selected genome's sequence dictionary and subset to the contigs the two share, producing the canonical `bams/{sample}.sorted.bam`. See [Reference contig filter](contig-filter.md).
