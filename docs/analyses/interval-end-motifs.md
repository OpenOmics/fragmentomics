# <code>finaletoolkit <b>interval-end-motifs</b></code>

## 1. About

> **Measures frequency of k-mer 5' end motifs in each region specified in a BED file and writes data into a table.**

The per-region counterpart of [`end-motifs`](end-motifs.md). Instead of one 256-entry frequency vector for the whole genome, it produces one vector per genomic window, so nuclease preference becomes a spatial signal. Regions where the end-motif profile departs from the genome average point to locally different chromatin accessibility or nuclease activity.

In this pipeline the windows are the fixed-size intervals tiling the genome, built at the start of each run at the width given by [`--interval`](../usage/run.md) (default `1mb`), so the output is a wide table: windows down the rows, 4-mer motifs across the columns.

## 2. As the pipeline runs it

```text
finaletoolkit interval-end-motifs bams/{sample}.sorted.bam <ref2bit> <intervals> \
    -q <--mapscore> \
    -k 4 \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -o interval_end_motifs/{sample}_endmotif_interval.tsv \
    -w <threads> \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `refseq_file`
> **2bit reference sequence.**
> *value:* the `ref2bit` file for the selected `--genome` build
>
> Motifs are read from the reference rather than from the reads, so base-calling errors at fragment ends do not affect the tally.

---
  `intervals`
> **Windows to tabulate motifs over.**
> *value:* `intervals/{genome}_{size}_intervals.bed`, tiled from the selected build's reference FastA at the [`--interval`](../usage/run.md) width

---
  `-k`
> **k-mer length.**
> *value:* `4` — 256 motifs, matching [`end-motifs`](end-motifs.md)

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `--mapscore` (default 20)

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)
>
> Passed as `-min`/`-max`. FinaleToolkit still accepts `-lo`/`-hi` as deprecated aliases for these; the pipeline uses the current names.

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `16`)

!!! note "Both strands are counted"

    `-B, --single-strand` would restrict the tally to one strand. The pipeline does not pass it, so both strands' 5' ends contribute — consistent with [`end-motifs`](end-motifs.md). Note the flag name differs between the two subcommands: `--single-strand` here, `--no-both-strands` there, for the same `-B` short form.

## 3. Output

`interval_end_motifs/{sample}_endmotif_interval.tsv` — a table with one row per interval and one column per 4-mer motif, holding that motif's frequency within the interval.

This is the largest per-sample text output in the pipeline: genome-wide windows against 256 motif columns. At the default `--interval 1mb` that is a few thousand rows; a narrow width such as `5kb` runs to hundreds of thousands. Its memory allocation is correspondingly the largest of the motif steps.

## 4. Controlling it

| What | How |
|------|-----|
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Which windows | `--interval`, which sets the width the windows are tiled at |
| Reference sequence | the `ref2bit` entry for the selected build |
| Threads / memory / walltime | the `interval_end_motifs` entry of `config/cluster.json` |

```json
"interval_end_motifs": {
    "threads": 16,
    "mem": "96G",
    "time": "08:00:00",
    "partition": "norm"
}
```

The 96 GB allocation reflects the size of the interval × motif table held in memory while it is assembled. If this step is killed for exceeding memory, either raise `mem` or pass a coarser `--interval` width.

## 5. Requires

`ref2bit` **and** `reference_fa` (which the interval BED is tiled from) for the selected genome build. Both bundled builds provide both.

## 6. Related

- [`end-motifs`](end-motifs.md) — the same measurement genome-wide
- [`mds`](mds.md) — motif diversity score, computed from the genome-wide vector
