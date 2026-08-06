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
3. **Coordinate sorts and converts to BAM** with `samtools sort`. CRAM and SAM inputs are converted to BAM as part of this step — there is no separate conversion step. Already coordinate-sorted input still passes through the sort, which is cheap relative to re-sorting and guarantees the ordering rather than trusting the header's `SO:` tag.
4. **Indexes the result**, producing a `.bai` alongside each BAM.
5. **Refreshes the index timestamp** so the index is never older than the BAM it indexes. Without this, `samtools` and `pysam` emit "index is older than data file" warnings on some filesystems where the two writes land in the same timestamp granularity.

The output lands in `staged_bams/`.

### 2.1 Duplicate sample names are rejected

Because sample names are derived from basenames, two inputs from different directories can collide:

```text
/path/one/Sample_A.bam   ─┐
                          ├─► both normalize to Sample_A
/path/two/Sample_A.cram  ─┘
```

`stage_bams` detects this before writing anything and fails with the list of colliding names. It does not silently pick one, and it does not disambiguate by adding a suffix — a collision means the two files would have produced results attributed to the same sample, and only the user knows which name each should have. Rename the inputs and re-run.

## 3. What is *not* done

`stage_bams` deliberately does not modify the read content of the input:

- **Duplicates are not marked or removed.** Unlike the [FastQ path](fastq-alignment.md), which marks and then drops duplicates, staged BAMs are taken as-is. If your BAMs have not been deduplicated, coverage will be inflated and fragment-length and end-motif distributions will be biased toward whatever the PCR amplified. Deduplicate before providing them.
- **No read filtering is applied.** The FastQ path keeps only properly-paired, primary, mapped, non-duplicate reads (`-f 2 -F 3852`). Staged BAMs retain every record they arrived with. FinaleToolkit applies its own per-analysis mapping-quality threshold (`-q 30`) and fragment-length window at read time, but it does not remove supplementary or secondary alignments for you.
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

`stage_bams` is a single step that loops over all inputs serially, so its walltime scales with the *total* size of your input set rather than the largest single file. For large cohorts, raise `time` in `config/cluster.json` in the output directory. The threads are passed to `samtools sort`.

## 5. Next step

Once staging completes, `staged_bams/{sample}.sorted.bam` is compared against the selected genome's sequence dictionary and subset to the contigs the two share, producing the canonical `bams/{sample}.sorted.bam`. See [Reference contig filter](contig-filter.md).
