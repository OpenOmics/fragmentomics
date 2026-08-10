# FastQ alignment

## 1. About

When the pipeline is given paired-end Illumina FastQ files, it runs the *full* pipeline: reads are aligned to the reference genome build selected with `--genome`, duplicates are marked, the alignment is filtered down to analysis-grade read pairs, and the result is written as the canonical `bams/{sample}.sorted.bam` that every downstream analysis consumes.

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

1. **Every downstream reference file is built against the same primary-only assembly.** The 2bit sequence, the 5 kb interval BED, the chrom.sizes file, the TSS BEDs, the gap track and the blacklist for a build all describe the same contig set. Aligning to the matching index means fragment coordinates, interval definitions and sequence lookups agree by construction.
2. **It is the reason BAM input needs a [contig filter](contig-filter.md) but FastQ input does not.** Reads aligned here can only land on contigs the reference declares, so a FastQ run produces a BAM whose contig set already matches the reference exactly. User-supplied BAMs carry whatever contig set their original aligner used, which is usually a full assembly with several hundred extra scaffolds.

Aligning to a build also selects the reference set used by the analysis steps. See [Reference files](references.md) for the full list of files each build provides.

## 3. What `align_fastq` does

The step is a single fused pipe — no intermediate BAM is written to disk between stages, which is what keeps the disk footprint and I/O of an alignment run manageable:

```text
${bwa_bin} mem -t <threads> -R "@RG\tID:<sample>\tSM:<sample>\tPL:ILLUMINA\tLB:<sample>" \
    <bwamem2_index> <sample>.R1.fastq.gz <sample>.R2.fastq.gz
  │
  ├─► samtools sort -n      query-name sort, required by fixmate
  ├─► samtools fixmate -m   fill in mate coordinates and add mate-score tags
  ├─► samtools sort         coordinate sort
  ├─► samtools markdup      flag optical/PCR duplicates
  └─► samtools view -f 2 -F 3852
         │
         ▼
    bams/<sample>.sorted.bam   (+ .bai)
```

### 3.1 Choosing the `bwa-mem2` binary

`bwa-mem2` is not distributed as one portable executable. The upstream release ships a separate binary per SIMD instruction set — `bwa-mem2.avx2`, `bwa-mem2.sse42`, `bwa-mem2.sse41` and others — each compiled for instructions that a CPU lacking them cannot execute at all. Running the wrong one does not degrade gracefully; it dies with `SIGILL` (illegal instruction).

So the binary is not hardcoded. The step's first action is to resolve `${bwa_bin}` by running `workflow/scripts/python_cpu_arch.py`, which reads the `flags` line of `/proc/cpuinfo` and prints the best-supported variant:

| Priority | CPU flag | Binary |
|---------|----------|--------|
| 1 | `avx2` | `bwa-mem2.avx2` |
| 2 | `sse4_2` | `bwa-mem2.sse42` |
| 3 | `sse4_1` | `bwa-mem2.sse41` |
| fallback | — | `bwa-mem2` (upstream dispatcher shim, decides for itself at runtime) |

Two properties of this are deliberate:

- **AVX-512 is never selected**, even on a CPU that advertises `avx512bw`. The upstream `bwa-mem2.avx512bw` binary has a long history of segfaults on both Intel and AMD parts, plus an AMD Zen–specific failure where the CPU reports `avx512bw` support but the binary takes Intel-tuned paths that fail at runtime (bwa-mem2 issues [#50](https://github.com/bwa-mem2/bwa-mem2/issues/50), [#112](https://github.com/bwa-mem2/bwa-mem2/issues/112), [#160](https://github.com/bwa-mem2/bwa-mem2/issues/160), [#199](https://github.com/bwa-mem2/bwa-mem2/issues/199), [#254](https://github.com/bwa-mem2/bwa-mem2/issues/254)). The marginal speedup is not worth a run that crashes hours in.
- **Detection happens on the compute node, not at submission.** It runs inside the step's shell block rather than when the workflow DAG is built, because in cluster mode the host that submits the job and the node that runs it need not share a CPU generation. Deciding at submit time on a newer login node would hand a partition's older nodes a binary they cannot run.

Detection failures are fatal before any alignment work starts: an unreadable `/proc/cpuinfo`, no CPU flags at all, or a non-x86 CPU (ARM — no `bwa-mem2` variant applies) each abort the step immediately, as does a resolved binary that is somehow absent from `PATH`. The chosen binary is echoed to the step's log, and it is also recorded in the output BAM's `@PG` header line, so after the fact you can confirm which variant produced a given alignment:

```bash
samtools view -H bams/<sample>.sorted.bam | grep '^@PG'
# @PG  ID:bwa-mem2  PN:bwa-mem2  VN:2.2.1  CL:bwa-mem2.avx2 mem -t 32 ...
```

Because the detection script runs inside the aligner container, that image provides `python3` alongside `bwa-mem2` and `samtools`.

### 3.2 Read groups

Each output BAM carries a single `@RG` line with `ID`, `SM` and `LB` all set to the sample name and `PL:ILLUMINA`. The sample name is taken from the FastQ basename, so per-sample provenance survives into the BAM header and into any downstream merge.

### 3.3 Duplicate marking

`samtools markdup` flags duplicates rather than removing them, and the subsequent filter drops them. The `-m` flag on `fixmate` is what makes this possible: it adds the mate-score tag `markdup` needs to pick which copy of a duplicate set to keep.

Because `markdup` runs inside a pipe without `-f`, it does not emit a metrics file. Duplicate rates are instead recovered from `samtools flagstat` by the [`bam_stats`](quality-control.md) step and surfaced in the MultiQC report.

### 3.4 Read filtering

Only reads passing both of the following survive into the analysis BAM:

  `-f 2`
> **Keep only properly-paired reads.**
>
> Flag `0x2`. The aligner marks a pair as "proper" when both mates align to the same contig in the expected orientation and within the expected distance. Fragment length, fragment midpoint, coverage and cleavage position are all inferred from the pair's coordinates, so an improperly-paired read has no meaningful fragment to contribute.

---
  `-F 3852`
> **Drop reads matching any of the following flags.**
>
> | Flag | Meaning | Why it is dropped |
> |------|---------|-------------------|
> | `0x4` (4) | read unmapped | no coordinates |
> | `0x8` (8) | mate unmapped | no fragment can be inferred |
> | `0x100` (256) | secondary alignment | would double-count a fragment |
> | `0x200` (512) | fails vendor QC | low-confidence base calls |
> | `0x400` (1024) | PCR/optical duplicate | inflates coverage and biases fragment-length and motif distributions |
> | `0x800` (2048) | supplementary alignment | chimeric segment, not an independent fragment |
>
> `4 + 8 + 256 + 512 + 1024 + 2048 = 3852`.

The net effect is that `bams/{sample}.sorted.bam` contains exactly one primary, non-duplicate, properly-paired alignment per surviving fragment — which is the input contract every FinaleToolkit analysis in this pipeline assumes. Note that no mapping-quality filter is applied here; each analysis applies its own `-q 30` threshold at read time, so the same BAM can be re-analyzed at a different threshold without re-aligning.

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

The `lscratch:200` request backs the temporary sort directories used by the two `samtools sort` stages; these are created under `--tmp-dir` and removed when the step exits, including on failure. Edit `config/cluster.json` in the output directory to change any of these for a given run.

## 5. Verifying the result

After a FastQ run, the [`bam_stats`](quality-control.md) step writes `samtools stats`, `flagstat` and `idxstats` reports for each analysis BAM into `qc/`, and MultiQC aggregates them into `multiqc/multiqc_report.html`. Useful checks there:

- **Total reads** in `flagstat` versus input FastQ read count — the gap is the reads removed by the filter above.
- **Duplicate rate** — high rates indicate a low-complexity library.
- **`idxstats` per-contig counts** — reads should be distributed across the primary chromosomes with no contig unexpectedly empty.
