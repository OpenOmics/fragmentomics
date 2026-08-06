# Reference contig filter

## 1. About

This step applies only to the **BAM input path**. It runs once per sample, after [BAM normalization](bam-normalization.md), and produces the canonical analysis BAM:

```text
staged_bams/{sample}.sorted.bam
         │
         │  filter_reference_contigs
         ▼
bams/{sample}.sorted.bam   +   qc/{sample}.contig_validation.txt
```

It exists because of a mismatch that is nearly universal for user-supplied alignments. This pipeline's references are **"clean" primary-assembly builds** — the primary chromosomes only. Real-world BAMs are usually aligned against a *full* assembly: the same primary chromosomes plus several hundred alt, decoy, patch and unplaced scaffolds. A BAM built against GRCh38 with decoys typically declares around 700 contigs where the reference declares 25.

The step does two things at once:

1. **Verifies** that the BAM really was aligned to the selected `--genome` build, by comparing contig names and lengths against the build's sequence dictionary.
2. **Filters** the BAM down to the contigs that the BAM and the reference share, writing a new BAM rather than failing on the extras.

!!! note

    The step is only defined when the selected genome build has a `dict` entry in `config/genome.json`. Both bundled builds (`hg19`, `hg38`) do. For a build without one there is nothing to compare against, so `stage_bams` writes `bams/` directly and this step is absent from the workflow — existing configurations keep working unchanged.

## 2. What is compared

The reference side of the comparison is a Picard-style sequence dictionary, whose `@SQ` records look like this:

```text
@SQ	SN:chr1	LN:248956422	M5:2648ae1bacce4ec4b6cf337dccf37e3b	UR:file:///data/OpenOmics/references/fragmentomics/hg38_clean.fa
```

**Only `SN` (contig name) and `LN` (contig length) are used.** `UR` and `M5` are deliberately ignored:

  `UR` — the URI of the FastA the dictionary was built from
> This is a path on the filesystem where the reference was *created*. It is essentially never present in an aligner-produced BAM header, and even when it is, it legitimately differs between the reference and the input for the same assembly — the two files live in different places. Comparing `UR` would report a mismatch for every correctly-aligned BAM, which makes it useless as a validation signal.

---
  `M5` — an MD5 checksum of the contig sequence
> `M5` is the *right* way to prove two contigs are byte-identical, but aligners do not populate it by default. In practice the field is absent from input BAM headers, so requiring it would fail every real input.

`SN` + `LN` is what remains reliably available on both sides, and it is sufficient to distinguish assemblies in practice: `chr1` is 248,956,422 bp in hg38 and 249,250,621 bp in hg19, so a build mismatch shows up immediately as a length disagreement.

The comparison is also **order-independent and count-independent**. The reference dictionary lists 25 contigs in lexical order; a full-assembly BAM lists ~700 in karyotypic order. Neither the ordering nor the total count is compared — only the per-name lengths of contigs that appear on both sides.

## 3. How each contig is classified

| Class | Condition | Action |
|-------|-----------|--------|
| **matched** | same `SN`, same `LN` in both | **kept** in the output BAM |
| **extra in BAM** | in the BAM, absent from the reference | **dropped** — this is the point of the step (alt/decoy/patch scaffolds) |
| **missing from BAM** | in the reference, absent from the BAM | **warned about**, output simply lacks the contig |
| **length mismatch** | same `SN`, different `LN` | **fatal** — the run stops |

### 3.1 Why a length mismatch is fatal

A shared contig name with a different length means the BAM was aligned to a *different assembly* than the one selected — an hg19 BAM run against the hg38 reference, for example. Every read coordinate in that BAM refers to a different position than the pipeline's interval, TSS and 2bit files assume. Subsetting contigs cannot repair that; the coordinates themselves are wrong. The step reports the first ten offending contigs with both lengths and exits non-zero:

```text
Fatal: contig length mismatch between the input BAM and the selected reference genome.
2 shared contig name(s) have different lengths, which means this BAM was aligned to a
different assembly:
    chr1: reference=248956422 bam=249250621
    chr2: reference=242193529 bam=243199373
Subsetting cannot correct this. Re-run with the matching --genome, or realign the input.
```

A BAM sharing **no** contigs at all with the reference is likewise fatal — usually a naming-convention difference (`1, 2, 3…` versus `chr1, chr2, chr3…`) rather than a wrong build, but in either case the input cannot be analyzed against this reference.

!!! warning "What this check cannot catch"

    A build difference that leaves all shared contig lengths identical is undetectable this way. The check proves the BAM is *consistent* with the selected build's contig geometry; it does not prove the sequences are identical. Only `M5` could do that, and aligners do not write it.

## 4. Why a new BAM instead of a header edit

Dropping contigs from a BAM is not a header operation, even though it looks like one. Two hazards make the naive version silently wrong, and the implementation handles both:

### 4.1 Reference indices, not names

**BAM records store their reference as an integer index into the `@SQ` list**, not as a contig name. Removing `@SQ` entries renumbers every remaining contig, so a record that pointed at index 20 now points at whatever contig ended up in slot 20. The result is a valid-looking BAM in which reads have been silently reassigned to the wrong chromosomes.

The filter avoids this by round-tripping the records through **SAM text**, which refers to contigs by *name*:

```text
reduced header (SAM text)  ─┐
                            ├─► samtools view -b ─► bams/{sample}.sorted.bam
headerless SAM records     ─┘
```

The final `samtools view -b` rebuilds the name-to-index mapping against the reduced header, so every read keeps the contig it started on.

### 4.2 Dangling mates

A read can sit on a kept contig while its mate sits on a dropped one. When the mate's contig disappears, `samtools` cannot resolve `RNEXT` and rewrites it to `*` — but it leaves `PNEXT` and the mate-mapped flag untouched. The result is a self-inconsistent record claiming a mapped mate at a position on no contig.

Such reads are excluded entirely, and the count is recorded in the report as `reads_dropped_dangling_mate`. This costs nothing analytically: a cross-contig pair is never a proper pair, and every FinaleToolkit analysis in this pipeline works from proper pairs, so these reads would have been ignored downstream regardless.

### 4.3 What is preserved

`@HD`, `@RG`, `@PG` and `@CO` header lines are carried through verbatim, so read groups and the provenance of the original alignment survive the subset. Kept `@SQ` records are emitted in **input-header order** rather than reference-dictionary order, which is what keeps a coordinate-sorted input coordinate-sorted on the way out. The output is indexed.

## 5. The validation report

Every run writes `qc/{sample}.contig_validation.txt`. It is a tab-delimited record of the comparison, written *before* the fatal checks so it survives as a diagnostic when validation fails:

```text
# Reference contig validation
bam	/data/$USER/output/staged_bams/Sample_A.sorted.bam
sequence_dictionary	/data/OpenOmics/references/fragmentomics/hg38_clean.dict
# Comparison uses @SQ SN and LN only; UR and M5 ignored.
reference_contigs	25
bam_contigs	706
matched	25
length_mismatches	0
missing_from_bam	0
extra_in_bam_dropped	681
reads_dropped_dangling_mate	0

# kept contigs (name, length)
kept	chr1	248956422
kept	chr2	242193529
...

# bam contigs absent from the reference (dropped)
dropped	chr1_KI270706v1_random
dropped	chrUn_GL000195v1
...
```

The report is an explicit target of the workflow, so it is produced even when every downstream output is already up to date. It is written into the shared `qc/` directory but is not parsed by MultiQC, which ignores it.

## 6. Resources

```json
"filter_reference_contigs": {
    "threads": 8,
    "mem": "16G",
    "time": "08:00:00",
    "gres": "lscratch:200",
    "partition": "norm"
}
```

The step streams rather than buffering, so memory is modest; walltime scales with the number of reads on kept contigs. It runs once per sample and samples run in parallel.

## 7. Reading the outcome

  **`extra_in_bam_dropped` is large (hundreds)**
> Expected and normal. Your BAM was aligned to a full assembly and the alt/decoy/patch scaffolds have been removed. Reads on those scaffolds are excluded from the analysis, which is the intended behavior — the pipeline's reference files describe only the primary assembly.

---
  **`missing_from_bam` is non-zero**
> The reference declares a contig your BAM's header does not. Often `chrM` or `chrY` (a female sample aligned against a Y-less reference, or a mitochondria-excluded alignment). The output simply lacks that contig, and any interval on it produces no data. Check the warning in the job log against your expectations for the sample.

---
  **`matched` is much smaller than `reference_contigs`**
> Your BAM covers only part of the assembly — a targeted panel, a single-chromosome test file, or a subset BAM. The analysis will run, but genome-wide features (particularly [DELFI](../analyses/delfi.md), which expects autosome-wide 100 kb bins) will be computed from a fraction of the genome.

---
  **The step failed with a length mismatch**
> Re-run with the `--genome` build your BAMs were actually aligned to, or realign them. Consult the `length_mismatch` lines in the report to identify which build you have: `chr1` at 249,250,621 bp is hg19; at 248,956,422 bp it is hg38.
