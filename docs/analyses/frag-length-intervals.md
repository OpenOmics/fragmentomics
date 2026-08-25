# <code>finaletoolkit <b>frag-length-intervals</b></code>

## 1. About

> **Retrieves fragment length summary statistics over intervals defined in a BED file.**

Where [`frag-length-bins`](frag-length-bins.md) gives one distribution for the whole genome, this gives summary statistics *per genomic window*: the mean, median, standard deviation, minimum and maximum fragment length in each interval. This turns fragment length from a single sample-level profile into a spatial signal, so regions where cfDNA is systematically shorter or longer than the genome average can be located.

In this pipeline the intervals are the fixed-size windows tiling the genome, built at the start of each run at the width given by [`--interval`](../usage/run.md) (default `1mb`).

## 2. As the pipeline runs it

```text
finaletoolkit frag-length-intervals bams/{sample}.sorted.bam <intervals> \
    -q <--mapscore> \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -p any \
    -o frag_length_intervals/{sample}_frag_interval.bed \
    -w <threads> \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `interval_file`
> **Windows to summarize over.**
> *value:* `intervals/{genome}_{size}_intervals.bed`, tiled from the selected `--genome` build's reference FastA at the [`--interval`](../usage/run.md) width
>
> See [Reference files](../pipeline/references.md).

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `--mapscore` (default 20)

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)
>
> Note that this truncates the distribution *before* the statistics are computed, so the reported mean and standard deviation describe fragments inside the window only.

---
  `-p, --intersect-policy`
> **How a fragment is assigned to an interval.**
> *value:* `any`
>
> A fragment contributes to every window it overlaps. A long fragment spanning a window boundary is therefore counted in both, which is the right behavior for characterizing the fragment population present at a locus. Double counting at boundaries matters more the narrower the `--interval` width is.

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `16`)

!!! note "Short-read threshold is left at its default"

    FinaleToolkit accepts `-s, --short-reads` to set the threshold defining the short-read fraction (default 150 bp). The pipeline does not override it. If you need short-to-long fragment ratios as a feature, [`delfi`](delfi.md) computes them explicitly and with GC correction.

## 3. Output

`frag_length_intervals/{sample}_frag_interval.bed` — a BED where each interval carries its fragment-length summary statistics: mean, median, standard deviation, minimum and maximum.

## 4. Controlling it

| What | How |
|------|-----|
| Whether this step runs | requires a non-zero `--split-interval` |
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Which windows | `--interval`, which sets the width the windows are tiled at |
| Threads / memory / walltime | the `frag_length_intervals` entry of `config/cluster.json` |

```json
"frag_length_intervals": {
    "threads": 16,
    "mem": "8G",
    "time": "08:00:00",
    "partition": "norm"
}
```

## 5. Requires

A non-zero `--split-interval` (default `5000`), and `reference_fa` for the selected build, which the interval BED is tiled from. Both bundled builds provide it.
