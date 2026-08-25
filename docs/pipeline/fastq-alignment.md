# FastQ alignment

## 1. About

When the pipeline is given paired-end Illumina FastQ files, it runs the *full* pipeline: reads are adapter-trimmed and quality-filtered, aligned to the reference genome build selected with `--genome`, the alignment is filtered down to analysis-grade read pairs — with any mate orphaned by that filter dropped alongside it — duplicates are marked (but kept), and the result is written as the canonical `bams/{sample}.sorted.bam` that every downstream analysis consumes. Read QC is collected on both sides of the alignment.

This is handled by a single workflow step, `align_fastq`, which runs once per sample.

!!! note

    Single-end FastQ data is **not** supported. All of the fragmentation features this pipeline computes are derived from inferred fragment coordinates, which require both mates of a properly-paired read.

## 2. What reference is aligned to

Reads are aligned with [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) against the pre-built index for the genome build chosen with `--genome`. The index path comes from the `bwamem2_index` entry of `config/genome.json`:

| `--genome` | `bwa-mem2` index |
|-----------|------------------|
| `hg38`    | `/data/OpenOmics/references/fragmentomics/hg38/bwamem2_index/hg38_clean.fa` |
| `hg19`    | `/data/OpenOmics/references/fragmentomics/hg19/bwamem2_index/hg19_clean.fa` |

Both are indexes of the **"clean" primary-assembly FastA** for that build — the primary chromosomes only, without alt, decoy or patch scaffolds. This matters for two reasons:

1. **Every downstream reference file is built against the same primary-only assembly.** The 2bit sequence, the chrom.sizes file, the TSS BEDs, the gap track and the blacklist for a build all describe the same contig set, and the genomic interval BED is tiled from that same FastA at run time. Aligning to the matching index means fragment coordinates, interval definitions and sequence lookups agree by construction.
2. **It is the reason BAM input needs a [contig filter](contig-filter.md) but FastQ input does not.** Reads aligned here can only land on contigs the reference declares, so a FastQ run produces a BAM whose contig set already matches the reference exactly. User-supplied BAMs carry whatever contig set their original aligner used, which is usually a full assembly with several hundred extra scaffolds.

Aligning to a build also selects the reference set used by the analysis steps. See [Reference files](references.md) for the full list of files each build provides.

## 3. What `align_fastq` does

FastQC runs first on the raw mates, `fastp` then trims adapters and drops low-quality pairs, the alignment itself is a single fused pipe — no intermediate BAM is written to disk between stages, which is what keeps the disk footprint and I/O of an alignment run manageable — and FastQC runs once more on the finished BAM:

```text
fastqc <sample>.R1.fastq.gz <sample>.R2.fastq.gz   pre-alignment read QC
  │
fastp --detect_adapter_for_pe \                    adapter trimming, then
      --average_qual <baseqscore>                  mean-base-quality filter
  │                                                (to $tmp, not kept)
${bwa_bin} mem -t <threads> -R "@RG\tID:<sample>\tSM:<sample>\tPL:ILLUMINA\tLB:<sample>" \
    <bwamem2_index> $tmp/<sample>.R1.trimmed.fastq.gz $tmp/<sample>.R2.trimmed.fastq.gz
  │                                                                   │
  ├─► samtools view -f 2 -F 2828 -q <mapscore>   flag + MAPQ filter   │
  ├─► samtools sort -n        query-name sort, what fixmate needs     │
  │                                          ── name-ordered ─────────┤
  ├─► samtools fixmate -m -r  fill mate coords, add ms tags,          │
  │                           de-pair mates left behind by the filter │
  ├─► samtools view -f 2      drop those de-paired orphans   ─────────┘
  ├─► samtools sort           coordinate sort
  └─► samtools markdup -f     flag optical/PCR duplicates, write dup report
         │
         ▼
    bams/<sample>.sorted.bam   (+ .bai)
         │
  fastqc <sample>.sorted.bam                       post-alignment read QC
```

Three properties of this ordering are deliberate:

- **Two sorts, each required by the stage it feeds.** `samtools fixmate` pairs records by walking adjacent reads of the same name, so it requires query-name-ordered input and `sort -n` is run to guarantee it rather than depending on the order `bwa-mem2` happens to emit. The coordinate sort that follows is equally non-optional: `markdup` and every downstream analysis need coordinate order. Both spill to the node's temporary directory, so the step's `lscratch` allocation covers two passes over the read data.
- **Filtering runs before the name sort,** so that mates orphaned by the filter can be cleaned up downstream and so the sort only handles reads that survive. See [§3.5](#35-orphaned-mates).
- **`markdup` runs last, after every filter.** See [§3.3](#33-duplicate-marking).

### 3.1 Choosing the `bwa-mem2` binary

`bwa-mem2` is not distributed as one portable executable. The upstream release ships a separate binary per SIMD instruction set — `bwa-mem2.avx2`, `bwa-mem2.sse42`, `bwa-mem2.sse41` and others — each compiled for instructions that a CPU lacking them cannot execute at all. Running the wrong one does not degrade gracefully; it dies with `SIGILL` (illegal instruction). The plain `bwa-mem2` binary is not a fourth variant but an upstream dispatcher shim that inspects the CPU at startup and re-execs the appropriate one.

So the binary is not hardcoded. The step's first action is to resolve `${bwa_bin}` by running `workflow/scripts/python_cpu_arch.py`, which reads the `vendor_id` line of `/proc/cpuinfo` and selects on **CPU vendor**:

| `vendor_id` | Binary | Why |
|-------------|--------|-----|
| `GenuineIntel` | `bwa-mem2` | Let the upstream dispatcher pick; its detection is tuned for Intel parts and gets AVX-512 right where it is safe. |
| `AuthenticAMD` | `bwa-mem2.avx2` | Pin AVX2 explicitly. Every AMD part this pipeline runs on supports it, and it side-steps the dispatcher's AVX-512 selection on Zen. |
| anything else | *(fatal)* | No `bwa-mem2` variant applies to a non-x86 CPU, so the step aborts rather than failing later with `SIGILL`. |

Two properties of this are deliberate:

- **AVX-512 is never selected on AMD.** The upstream `bwa-mem2.avx512bw` binary has a long history of segfaults on both Intel and AMD parts, plus an AMD Zen–specific failure where the CPU reports `avx512bw` support but the binary takes Intel-tuned paths that fail at runtime (bwa-mem2 issues [#50](https://github.com/bwa-mem2/bwa-mem2/issues/50), [#112](https://github.com/bwa-mem2/bwa-mem2/issues/112), [#160](https://github.com/bwa-mem2/bwa-mem2/issues/160), [#199](https://github.com/bwa-mem2/bwa-mem2/issues/199), [#254](https://github.com/bwa-mem2/bwa-mem2/issues/254)). Pinning AVX2 on AMD avoids that path entirely; the marginal speedup is not worth a run that crashes hours in.
- **Detection happens on the compute node, not at submission.** It runs inside the step's shell block rather than when the workflow DAG is built, because in cluster mode the host that submits the job and the node that runs it need not share a CPU generation or vendor. Deciding at submit time on a login node would hand a partition's nodes a binary chosen for the wrong hardware.

Detection failures are fatal before any alignment work starts: a missing `/proc/cpuinfo`, or a `vendor_id` that is neither Intel nor AMD (ARM, for instance), aborts the step immediately with an explanatory message, as does a resolved binary that is somehow absent from `PATH`. The chosen binary is echoed to the step's log, and it is also recorded in the output BAM's `@PG` header line, so after the fact you can confirm which variant produced a given alignment:

```bash
samtools view -H bams/<sample>.sorted.bam | grep '^@PG'
# @PG  ID:bwa-mem2  PN:bwa-mem2  VN:2.2.1  CL:bwa-mem2.avx2 mem -t 32 ...
```

Because the detection script runs inside the aligner container, that image provides `python3` alongside `bwa-mem2` and `samtools`.

### 3.2 Read groups

Each output BAM carries a single `@RG` line with `ID`, `SM` and `LB` all set to the sample name and `PL:ILLUMINA`. The sample name is taken from the FastQ basename, so per-sample provenance survives into the BAM header and into any downstream merge.

### 3.3 Duplicate marking

`samtools markdup` flags duplicates rather than removing them, and it is the **last** stage of the pipe — downstream of both `samtools view` filters. The ordering is what preserves the flags: a filter that ran after marking could act on `0x400`, and the conventional `-F 3852` would do exactly that. Running `markdup` last means nothing in this step can drop a read for being a duplicate, so the analysis BAM contains every properly-paired primary read with duplicates *identified* rather than deleted, and each downstream analysis is free to include or exclude them.

The `-m` flag on `fixmate` is what makes marking possible at all: it adds the mate-score tag `markdup` needs to pick which copy of a duplicate set to keep. Since `fixmate` now runs after the quality filter ([§3.5](#35-orphaned-mates)), those scores describe only the reads that survive it.

`markdup` runs with `-f`, so it writes a duplicate report to `qc/{sample}.markdup.stats.txt` — read/written/excluded/examined counts, paired and single duplicate counts, optical duplicates, and an estimated library size. MultiQC parses it natively, so it lands in the aggregate report alongside the `samtools` reports; the duplicate rate is also recoverable independently from `samtools flagstat` via [`bam_stats`](quality-control.md).

!!! warning "Duplicates reach the analysis steps"

    Because duplicates are no longer removed here, they are present in `bams/{sample}.sorted.bam` and count toward coverage, fragment-length and motif distributions unless the analysis step itself excludes `0x400`. A high duplicate rate in the markdup report is therefore worth acting on rather than just noting.

### 3.3.1 Read QC (`fastqc`)

FastQC runs twice per sample, both times writing into `qc/`:

| When | Input | Reports |
|------|-------|---------|
| Before alignment | `inputs/{sample}.R1.fastq.gz`, `inputs/{sample}.R2.fastq.gz` | `qc/{sample}.R1_fastqc.{zip,html}`, `qc/{sample}.R2_fastqc.{zip,html}` |
| After alignment | `bams/{sample}.sorted.bam` | `qc/{sample}.sorted_fastqc.{zip,html}` |

The report names are derived by FastQC from its input filenames, not chosen by the rule. The pre-alignment pair runs on the raw mates, *upstream* of `fastp`, so it describes the library exactly as sequenced — untrimmed, with its adapter content intact. The post-alignment report describes the reads that actually survived trimming, filtering and alignment, so comparing the two shows what the step removed. Both are parsed by MultiQC — the pre-alignment reports appear there as `{sample}.R1` and `{sample}.R2`, the post-alignment one as `{sample}`.

### 3.4 Adapter trimming and read filtering

`fastp` runs once, before alignment, and does two things: it trims adapter read-through out of the reads, and it drops pairs whose base quality is too low. `samtools view` then drops alignments that are not analysis-grade. Both quality thresholds come from the command line — [`--baseqscore`](../usage/run.md#22-analysis-options) and [`--mapscore`](../usage/run.md#22-analysis-options) — and both are inclusive lower bounds, so a read is kept when its score is greater than or equal to the value given. Passing `0` for either disables that filter.

#### Adapter trimming and base quality, before alignment (`fastp`)

```text
fastp --detect_adapter_for_pe \
      --average_qual <baseqscore> --unqualified_percent_limit 100 ...
```

**Adapter trimming.** cfDNA fragments are short — the mononucleosome peak sits near 167 bp — so a substantial fraction of any 2×150 bp library sequences straight off the end of the insert and into the sequencing adapter. `fastp` removes that read-through. For paired-end input its primary mechanism is **overlap analysis** rather than matching a known adapter sequence: the two mates of a short fragment overlap, the overlap reveals where the insert ends, and everything past it is trimmed. This suits the case that matters here, because a read only contains adapter when the fragment is shorter than the read length, and such a pair overlaps almost completely. `--detect_adapter_for_pe` additionally turns on auto-detection of the adapter sequence itself, which covers pairs that overlap too little for the overlap method to resolve; `fastp` disables it for paired-end input by default.

Trimming before alignment rather than relying on `bwa-mem2` to soft-clip means the fragment coordinates every downstream analysis is built on are inferred from insert sequence only.

**Quality filtering.** A read pair is discarded when its mean base quality is below `--baseqscore` (default 20). Filtering here rather than after alignment means failing pairs never consume alignment time at all.

Order matters between the two, and it is `fastp`'s: trimming happens first, so `--average_qual` is measured on the **trimmed** read. Adapter bases — which are frequently the lowest-quality bases in the read, being at the 3′ end — therefore cannot drag a good pair below the threshold. Verified against `fastp 0.23.4`: pairs of 150 bp reads carrying 50 bp of adapter read-through at Q2 over 100 bp of insert at Q40 have an untrimmed mean of 27.3 and a trimmed mean of 40, and they pass `--average_qual 35`.

The trimmed mates are written to the step's temporary directory and consumed directly by the aligner, so no preprocessed FastQ survives the step; `fastp`'s own before/after report is kept at `qc/{sample}.fastp.json` (plus an HTML copy) and is parsed by MultiQC, which is where the trimmed-read and adapter counts are visible.

Two `fastp` behaviors are worth calling out explicitly:

- `--unqualified_percent_limit 100` — overrides a default, disabling `fastp`'s separate per-base rule (by default, discarding a read when more than 40% of its bases fall below a per-base cutoff), so mean base quality is the only *quality* criterion applied and `--baseqscore` means exactly one thing. That is what makes `--baseqscore` mean the same thing on both input paths: the [BAM path](bam-normalization.md) enforces the identical mean-quality threshold with `samtools`, though it cannot trim.
- `--length_required 15` — left at its default, and only reachable now that trimming is on: a pair is dropped if either read trims shorter than 15 bp. Such a read is an adapter dimer or near-dimer with no usable insert. Note that this and the quality filter both act on the **pair**, so neither can orphan a mate.

#### Alignment flags and mapping quality (`samtools view`)

Only reads passing all of the following survive into the analysis BAM. The filter runs *before* duplicate marking, so `0x400` is deliberately absent from the exclusion mask, and *before* `fixmate`, so that the mates it strands can be dealt with ([§3.5](#35-orphaned-mates)):

  `-f 2`
> **Keep only properly-paired reads.**
>
> Flag `0x2`. The aligner marks a pair as "proper" when both mates align to the same contig in the expected orientation and within the expected distance. Fragment length, fragment midpoint, coverage and cleavage position are all inferred from the pair's coordinates, so an improperly-paired read has no meaningful fragment to contribute.

---
  `-F 2828`
> **Drop reads matching any of the following flags.**
>
> | Flag | Meaning | Why it is dropped |
> |------|---------|-------------------|
> | `0x4` (4) | read unmapped | no coordinates |
> | `0x8` (8) | mate unmapped | no fragment can be inferred |
> | `0x100` (256) | secondary alignment | would double-count a fragment |
> | `0x200` (512) | fails vendor QC | low-confidence base calls |
> | `0x800` (2048) | supplementary alignment | chimeric segment, not an independent fragment |
>
> `4 + 8 + 256 + 512 + 2048 = 2828`.
>
> This is the conventional `3852` **minus `0x400`** (PCR/optical duplicate). Duplicates are not marked yet at this point in the pipe, and by design they are never dropped — see [§3.3](#33-duplicate-marking).

---
  `-q <mapscore>`
> **Drop alignments below the mapping-quality threshold.**
>
> From `--mapscore`, default 20; a read is kept when its `MAPQ >= mapscore`. Multi-mapping and ambiguously-placed reads cannot be assigned a trustworthy fragment coordinate, which is what every fragmentation feature is derived from.
>
> The same value is passed to every FinaleToolkit analysis as its `-q` threshold, so one setting governs both what is written into the analysis BAM and what the analyses read out of it. See [shared conventions](../analyses/index.md#3-shared-conventions).

The net effect is that `bams/{sample}.sorted.bam` contains one primary, properly-paired alignment per surviving fragment, above the requested mapping-quality threshold, with PCR/optical duplicates flagged but retained, and with both mates of every pair present.

### 3.5 Orphaned mates

The filter above acts on individual records, not on pairs. A pair whose R1 has `MAPQ 60` and whose R2 has `MAPQ 5` loses only R2 at `--mapscore 20`, and the surviving R1 is then an **orphan**: it still carries `0x2` (properly paired) and mate coordinates pointing at a record that is no longer in the file. Every fragmentation feature this pipeline computes is derived from a pair's coordinates, so an orphan is at best dead weight and at worst a fragment inferred from a mate that does not exist.

Three stages cooperate to remove them, which is the reason the filter runs where it does:

```text
samtools view -f 2 -F 2828 -q <mapscore>   ── strands one mate of some pairs
samtools sort -n                           ── name order, so mates are adjacent
samtools fixmate -m -r                     ── sees the survivor as a singleton,
                                               clears 0x1/0x2 and the mate fields
samtools view -f 2                          ── drops what fixmate de-paired
```

`fixmate` can only do this because the stream reaches it in query-name order and before the coordinate sort — mates are adjacent, so a missing partner is visible. Once the BAM is coordinate-sorted the two mates are megabases apart and `fixmate` cannot pair them at all; that is why this guard exists on the FastQ path but not on the [BAM path](bam-normalization.md#22-read-filtering), where the input arrives already coordinate-sorted and re-collating it would mean an extra whole-file pass in a step whose only job is staging.

The three stages have distinct jobs, and the last one is what makes the guarantee:

  `fixmate -m`
> **Adds the `ms` mate-score tag** that `markdup` needs to choose which copy of a duplicate set to keep. Because it now runs *after* the filter, those scores are computed over the reads that actually survive, not over reads that are about to be discarded.

---
  `fixmate -r`
> **Removes unmapped and secondary leftovers.** These are already excluded by `-F 2828`, so in this pipe it is a no-op in the common case; it also removes some orphans outright, but not deterministically — a singleton in the middle of the stream is de-paired rather than deleted.

---
  `samtools view -f 2`
> **Drops the de-paired records.** This is the stage that actually removes orphans. It re-checks `0x2` only — no `-q`, no `-F` — because everything else was already applied upstream and re-applying it would be wasted work.

Verified against `samtools 1.13`: a pair whose second mate is filtered out has its survivor rewritten from flag `99` to flag `64` by `fixmate` (paired, proper-pair and mate-reverse bits cleared, `RNEXT`/`PNEXT`/`TLEN` zeroed), and flag `64` then fails `-f 2` and is dropped.

## 4. Resources

`align_fastq` is by far the most expensive step in the pipeline. Its defaults in `config/cluster.json`:

```json
"align_fastq": {
    "threads": 32,
    "mem": "200000m",
    "time": "72:30:00",
    "gres": "lscratch:200",
    "partition": "norm"
}
```

The `lscratch:200` request backs the temporary directories used by the two `samtools sort` stages (name and coordinate) and by `markdup` and FastQC, and it also holds the adapter-trimmed FastQ pair that `fastp` hands to the aligner (roughly the size of the input FastQ files); these are created under `--tmp-dir` and removed when the step exits, including on failure. Edit `config/cluster.json` in the output directory to change any of these for a given run.

## 5. Verifying the result

After a FastQ run, `qc/` holds this step's FastQC, `fastp` and markdup reports plus the `samtools stats`, `flagstat` and `idxstats` reports written by [`bam_stats`](quality-control.md), and MultiQC aggregates all of them into `multiqc/multiqc_report.html`. Useful checks there:

- **Total reads** in `flagstat` versus input FastQ read count — the gap is the pairs `fastp` dropped on base quality or post-trim length plus the reads removed by the filter above.
- **Singletons** in `flagstat` should be `0`, and *properly paired* should equal *total*. Anything else means an orphan reached the BAM, which the `fixmate` guard in [§3.5](#35-orphaned-mates) is there to prevent.
- **Duplicate rate** — from `qc/{sample}.markdup.stats.txt` (or `flagstat`); high rates indicate a low-complexity library, and note that those duplicates remain in the analysis BAM.
- **Pre- versus post-alignment FastQC** — per-base quality and adapter content on the raw mates, against the same metrics on the reads that survived trimming and filtering. Adapter content should be present in the pre-alignment reports and gone from the post-alignment one; FastQC runs on the *raw* mates, not on `fastp`'s output, so the "before" picture is genuinely untrimmed.
- **Adapter and trimming counts** in the `fastp` section of MultiQC — how much read-through was removed. A library showing no adapter at all in the pre-alignment FastQC and none trimmed here is either already trimmed or has an insert size comfortably above the read length.
- **`idxstats` per-contig counts** — reads should be distributed across the primary chromosomes with no contig unexpectedly empty.
