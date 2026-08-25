# BAM normalization

## 1. About

When the pipeline is given ready-made alignments instead of FastQ files, it runs the *full pipeline minus alignment*. The provided files still need to be brought into a uniform shape before the analysis steps can read them, because FinaleToolkit requires a coordinate-sorted, indexed BAM and user-supplied alignments arrive in whatever form their producer left them in.

Normalization is a three-stage process:

```text
<user-provided>.bam
         │
         │  stage_bams                 ← this page, §2
         ▼
staged_bams/{sample}.sorted.bam
         │
         │  remove_orphan_reads        ← this page, §3
         ▼
staged_bams/{sample}.paired.bam
         │
         │  filter_reference_contigs   ← see "Reference contig filter"
         ▼
bams/{sample}.sorted.bam
```

This page covers the first two stages. The third is documented in [Reference contig filter](contig-filter.md).

## 2. What `stage_bams` does

`stage_bams` runs once per pipeline invocation and processes every input file. For each one it:

1. **Determines the format from the file extension** — `.bam`, case-insensitively. A file whose extension is none of these is a hard error; the pipeline does not attempt to sniff the format from the file contents.
2. **Filters reads on mapping and base quality** with `samtools view`, using the thresholds given on the command line. See [§2.2](#22-read-filtering).
3. **Coordinate sorts and converts to BAM** with `samtools sort`. Already coordinate-sorted input still passes through the sort, which is cheap relative to re-sorting and guarantees the ordering rather than trusting the header's `SO:` tag. The filter above is piped straight into the sort, so no intermediate alignment file is written to disk.
4. **Indexes the result**, producing a `.bai` alongside each BAM.
5. **Refreshes the index timestamp** so the index is never older than the BAM it indexes. Without this, `samtools` and `pysam` emit "index is older than data file" warnings on some filesystems where the two writes land in the same timestamp granularity.

The output lands in `staged_bams/`.

### 2.1 Duplicate sample names are rejected

Because sample names are derived from basenames, two inputs from different directories can collide:

```text
/path/one/Sample_A.bam   ─┐
                          ├─► both normalize to Sample_A
/path/two/Sample_A.bam   ─┘
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

Both filters act on individual records, not on pairs, so either can keep one mate and drop the other. That is what the next step exists to repair — see [§3](#3-what-remove_orphan_reads-does).

### 2.3 Resources

```json
"stage_bams": {
    "threads": 48,
    "mem": "96G",
    "time": "12:00:00",
    "partition": "norm"
}
```

`stage_bams` is a single step that loops over all inputs serially, so its walltime scales with the *total* size of your input set rather than the largest single file. For large cohorts, raise `time` in `config/cluster.json` in the output directory. The threads are passed to both `samtools view` and `samtools sort`, which run concurrently as two ends of a pipe.

## 3. What `remove_orphan_reads` does

`stage_bams` filters record by record, so a pair whose R1 has `MAPQ 60` and whose R2 has `MAPQ 5` loses only R2 at `--mapscore 20`. The surviving R1 is then an **orphan**: it still carries `0x2` (properly paired) and mate coordinates pointing at a record that is no longer in the file. Nothing in the staged BAM marks it as such — `flagstat` still counts it as properly paired, because the flags on the record itself are unchanged. Reads whose mate went unmapped are the same problem arriving from the input BAM rather than being created by the filter, and `flagstat` does report those, as **singletons**.

Every fragmentation feature in this pipeline is derived from a pair's coordinates, so an orphan is at best dead weight and at worst a fragment inferred from a mate that does not exist. `remove_orphan_reads` removes both classes, one sample at a time, with the same pipe shape the [FastQ path](fastq-alignment.md#35-orphaned-mates) uses for the same purpose:

```text
samtools view -b -f 1 -F 2316   keep primary paired records, both mates mapped
  │
samtools sort -n                query-name order, so mates become adjacent
  │
samtools fixmate                recompute mate coordinates/ISIZE; de-pair any
  │                             read whose mate is missing (clears 0x1/0x2)
samtools view -b -f 1           drop those de-paired records
  │
samtools sort                   back to coordinate order
  │
  ▼
staged_bams/{sample}.paired.bam   (+ .bai, both deleted once the contig
                                   filter has consumed them)
```

Three details are worth knowing:

- **The `-F 2316` mask excludes unmapped (`0x4`), mate-unmapped (`0x8`), secondary (`0x100`) and supplementary (`0x800`) records.** The first two are the singletons — both members of such a pair go, the unmapped one on `0x4` and its partner on `0x8`. The other two are excluded because they are what makes the pairing question ill-defined: a secondary or supplementary record shares its read name with a primary one, so "every name appears exactly twice" cannot hold while they are present, and `fixmate` has no way to tell which record is the mate. This is the one place where the BAM path now removes records on flags; see [§4](#4-what-is-not-done).
- **`0x2` and `0x400` are deliberately *not* in the mask.** Discordant pairs — both mates present, but not flagged as a proper pair — are kept, since dropping them is a separate editorial decision and not one this step is making. Duplicate flags are left exactly as the input BAM carried them, so nothing about duplicate marking changes here.
- **`fixmate` is run without `-m`.** That flag adds the mate-score tag `samtools markdup` needs, and this path does not mark duplicates.

Single-end input is passed through unchanged rather than emptied. A BAM with no paired records has no orphans to remove but would not survive `view -f 1`, so the step checks the read count of its own output and copies the staged BAM through with a warning if the pipe removed everything. Both counts are read from the BAM indices rather than from a pass over the reads, so the check costs nothing on a normal paired-end run.

Verified against `samtools 1.13` on a staged BAM carrying 710 orphans among 555,122 records: 554,412 records survive, `flagstat` reports 0 singletons, and every surviving read name appears exactly twice.

### 3.1 Resources

```json
"remove_orphan_reads": {
    "threads": 16,
    "mem": "32G",
    "time": "12:00:00",
    "gres": "lscratch:200",
    "partition": "norm"
}
```

Unlike `stage_bams`, this step is one job per sample, so samples are processed in parallel. It sorts the staged BAM twice — once by name and once by coordinate — and both spill into `--tmp-dir`, which is what the `lscratch` request is for. Raise it for unusually deep BAMs; the temporary files are on the order of the BAM's own size and are removed when the step exits, including on failure.

## 4. What is *not* done

Apart from the quality filtering in [§2.2](#22-read-filtering) and the pairing cleanup in [§3](#3-what-remove_orphan_reads-does), the BAM path deliberately does not modify the read content of the input:

- **Duplicates are not marked.** The [FastQ path](fastq-alignment.md) runs `samtools markdup` so duplicates are at least *flagged* (it does not remove them either); staged BAMs are taken as-is, with whatever `0x400` flags — or none — they arrived with. If your BAMs have not been deduplicated, coverage will be inflated and fragment-length and end-motif distributions will be biased toward whatever the PCR amplified. Deduplicate before providing them.
- **Proper pairing is not required.** The FastQ path keeps only *properly-paired* reads (`-f 2 -F 2828`); this path requires a complete, mapped, primary pair (`-f 1 -F 2316`) but not the `0x2` flag, so discordant pairs reach the analysis BAM. FinaleToolkit re-applies the mapping-quality threshold and its fragment-length window at read time, which is what excludes the implausible fragment lengths such a pair implies.
- **Read groups are not added or rewritten.** Whatever `@RG` lines the input carries are preserved. So are `@PG` lines, which means the provenance of the original alignment survives into the analysis BAM.
- **Reads are not re-aligned.** The alignment coordinates in a staged BAM are exactly the ones its original aligner produced, against whatever reference that aligner used. The pipeline cannot detect a *wrong* reference from the coordinates alone — that is what the [contig filter](contig-filter.md) is for.

!!! warning "Staged BAMs are assumed to be aligned to the selected `--genome` build"

    `--genome` is required for BAM input, because it selects the interval, TSS, chrom.sizes, 2bit and blacklist files used by the analysis steps. If your BAMs were aligned to a different build, those reference files describe different coordinates than your reads do and the results will be meaningless. The [contig filter](contig-filter.md) catches the common, detectable version of this mistake — a contig present in both with a different length — and fails the run. It cannot catch a build difference that leaves all shared contig lengths identical.

## 5. Next step

Once the pairing cleanup completes, `staged_bams/{sample}.paired.bam` is compared against the selected genome's sequence dictionary and subset to the contigs the two share, producing the canonical `bams/{sample}.sorted.bam`. See [Reference contig filter](contig-filter.md).

When the selected build ships no sequence dictionary there is nothing to validate against, so that step does not exist and `remove_orphan_reads` writes `bams/{sample}.sorted.bam` itself. Both bundled builds, `hg19` and `hg38`, do ship one.
