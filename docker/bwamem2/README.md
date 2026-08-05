# bwa-mem2 aligner image

Alignment image for the fragmentomics **FastQ input path**. When `--input`
resolves to paired-end Illumina FastQ files, the `align_fastq` rule uses this
image to align reads to the reference genome and produce the same
`bams/{sample}.sorted.bam` that the BAM input path stages directly.

Contents:
- `bwa-mem2` v2.2.1 (paired-end alignment)
- `samtools` (fixmate, sort, markdup, proper-pair/mapped filtering, index)

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
docker build --platform linux/amd64 -t rroutsong/fragmentomics_bwamem2:0.0.1 .
docker push rroutsong/fragmentomics_bwamem2:0.0.1
```

The image URI is registered in `config/containers.json` under the `bwamem2`
key and can be cached locally with `fragmentomics cache`.
