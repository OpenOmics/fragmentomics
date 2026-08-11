# bwa-mem2 aligner image

Alignment image for the fragmentomics **FastQ input path**. When `--input`
resolves to paired-end Illumina FastQ files, the `align_fastq` rule uses this
image to align reads to the reference genome and produce the same
`bams/{sample}.sorted.bam` that the BAM input path stages directly.

It is also used by `bam_to_bed`, which runs on **both** input paths: that rule
needs `bedtools`, `bgzip` and `tabix`, and this is the image that carries them.

Contents:
- `bwa-mem2` v2.2.1 (paired-end alignment), including the per-instruction-set
  binaries the upstream release ships: `bwa-mem2.avx2`, `bwa-mem2.sse42`,
  `bwa-mem2.sse41`, `bwa-mem2.avx`, `bwa-mem2.avx512bw`, and the plain
  `bwa-mem2` dispatcher shim
- `samtools` (fixmate, sort, markdup, proper-pair/mapped filtering, index)
- `fastqc` — read QC, run by `align_fastq` on the raw mates before alignment
  and on the analysis BAM after it
- `fastp` — read preprocessing (trimming/filtering)
- `bedtools` — interval utilities, including the `bamtobed` conversion behind
  `bam_to_bed`
- `tabix` and `bgzip` — BGZF compression and interval indexing for
  `bam_to_bed`. Both come from Ubuntu's `tabix` package, *not* from `samtools`,
  so the package has to stay in the install list on its own for `beds/` output
  to be produced at all
- `python3` — required by `workflow/scripts/python_cpu_arch.py`, which
  `align_fastq` runs inside this image to pick the SIMD variant matching the
  compute node's CPU (see
  [FastQ alignment §3.1](../../docs/pipeline/fastq-alignment.md)). Do not drop
  it from the image: without an interpreter the rule cannot resolve a binary
  and every FastQ run fails at its first step.

## Reference index

`bwa-mem2` requires an index built from the reference FASTA. The index prefix
is configured per genome build in `config/genome.json` as `bwamem2_index`, and
the un-indexed FASTA as `reference_fa`. Both supported builds are self-hosted
under the fragmentomics reference directory, each in its own `bwamem2_index/`
subdirectory so the index components never collide with the finaletoolkit
references:

```
hg38  reference_fa   /data/OpenOmics/references/fragmentomics/hg38_clean.fa
      bwamem2_index  /data/OpenOmics/references/fragmentomics/hg38/bwamem2_index/hg38_clean.fa
hg19  reference_fa   /data/OpenOmics/references/fragmentomics/hg19_clean.fa
      bwamem2_index  /data/OpenOmics/references/fragmentomics/hg19/bwamem2_index/hg19_clean.fa
```

Each index expands the prefix to `<prefix>.{0123,amb,ann,bwt.2bit.64,pac}`. The
`_clean.fa` primary contigs (chr1–22, X, Y, M) match that build's finaletoolkit
references (`{build}.2bit`, `{build}.chrom.sizes`, and the TSS/interval BEDs)
exactly, so aligned BAM coordinates line up with every downstream analysis step.

### (Re)building an index

Indexes are built with the pinned bwa-mem2 version via SLURM (the build peaks
at a large memory footprint, so do not run it on a memory-capped shell). The
submission script lives next to the references:

```bash
sbatch /data/OpenOmics/references/fragmentomics/build_bwamem2_indexes.sbatch
```

To (re)build a single FASTA by hand on a suitably-resourced node:

```bash
module load bwa-mem2/2.2.1   # or run inside this image
bwa-mem2 index \
    -p /data/OpenOmics/references/fragmentomics/hg38/bwamem2_index/hg38_clean.fa \
    /data/OpenOmics/references/fragmentomics/hg38_clean.fa
```

## Build & push

```bash
docker build --platform linux/amd64 -t rroutsong/fragmentomics_bwamem2:0.0.4 .
docker push rroutsong/fragmentomics_bwamem2:0.0.4
```

The image URI is registered in `config/containers.json` under the `bwamem2`
key and can be cached locally with `fragmentomics cache`.
