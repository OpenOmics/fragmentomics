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
is configured per genome build in `config/genome.json` as `bwamem2_index`
(equal to the `reference_fa` path). Build it once, next to the FASTA, with the
**same** bwa-mem2 version pinned in the Dockerfile:

```bash
module load bwa-mem2/2.2.1   # or run inside this image
bwa-mem2 index /data/OpenOmics/references/fragmentomics/hg38_clean.fa
bwa-mem2 index /data/OpenOmics/references/fragmentomics/hg19_clean.fa
```

This writes `<prefix>.{0123,amb,ann,bwt.2bit.64,pac}` alongside the FASTA.

## Build & push

```bash
docker build --platform linux/amd64 -t rroutsong/fragmentomics_bwamem2:0.0.1 .
docker push rroutsong/fragmentomics_bwamem2:0.0.1
```

The image URI is registered in `config/containers.json` under the `bwamem2`
key and can be cached locally with `fragmentomics cache`.
