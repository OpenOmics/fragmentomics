# <code>finaletoolkit <b>wps</b></code>

## 1. About

> **Calculates Windowed Protection Score (WPS) over intervals defined in a BED file.**

WPS, from Snyder *et al.* (2016), measures nucleosome protection at base-pair resolution. For each position, it counts the fragments that **completely span** a window centered on it, minus the fragments whose **endpoint falls inside** that window:

```text
WPS(position) = (fragments spanning the window) − (fragments ending within the window)
```

A position protected by a nucleosome accumulates spanning fragments and few endpoints, giving a high score. A position exposed to nuclease attack accumulates endpoints, giving a low score. Run across a transcription start site, the resulting track reveals nucleosome positioning — and because nucleosome positioning is cell-type-specific, WPS around TSSs carries information about which tissues contributed the cfDNA.

This pipeline computes **L-WPS** (long-fragment WPS), the nucleosome-scale variant.

## 2. As the pipeline runs it

```text
finaletoolkit wps bams/{sample}.sorted.bam <tss> \
    -i <--split-interval> \
    -W 120 \
    -min 120 \
    -max 180 \
    -q 30 \
    -o wps/{sample}_wps_out_tss.bw \
    -w <threads> \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `site_bed`
> **Sites to compute WPS around.**
> *value:* the `tss` file for the selected `--genome` build
>
> Individual TSS positions. FinaleToolkit requires this file to be **sorted by contig then start**; the bundled `tss.hg38_sorted.bed` and `tss.hg19_sorted.bed` are. An unsorted site file produces incorrect output rather than an error.

---
  `-i, --interval-size`
> **Size of the window built around each site.**
> *value:* `--split-interval` (default `5000`)
>
> Each TSS is expanded to an interval of this size, centered on the site, and WPS is computed across every base of it. At the default this gives ±2.5 kb around each TSS — enough to capture the nucleosome-depleted region and several flanking nucleosomes.

---
  `-W, --window-size`
> **Sliding window used to compute each score.**
> *value:* `120`
>
> This is the window a fragment must span to count as protective. 120 bp approximates the nucleosome footprint, which is what makes the score read out nucleosome occupancy. Distinct from `-i`, which is the *extent* of the region scored.

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `120` / `180` — **fixed, not taken from `--fragment-minimum`/`--fragment-maximum`**
>
> This restriction to 120–180 bp is what defines **L-WPS**: only mononucleosome-sized fragments contribute. Including shorter fragments would mix in a sub-nucleosomal signal that reflects a different protection mechanism and would blur the nucleosome periodicity. These are FinaleToolkit's own defaults for the subcommand, and the pipeline states them explicitly rather than relying on them.
>
> !!! note
>
>     This is the one analysis in the pipeline that ignores the `--fragment-minimum` and `--fragment-maximum` options. Widening those has no effect on WPS.

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `30`

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `24`)

!!! note "`--chrom-sizes` is not passed here"

    `wps` accepts an optional `-c, --chrom-sizes`; the pipeline does not pass it, since the site BED bounds the regions scored. The chrom.sizes file is required by [`adjust-wps`](adjust-wps.md) downstream, which does receive it.

## 3. Output

`wps/{sample}_wps_out_tss.bw` — a bigWig of the raw per-base WPS across every TSS window.

Raw WPS is noisy and carries a long-wavelength trend that obscures the nucleosome periodicity, so it is normally not interpreted directly. Two things consume it:

- [`adjust-wps`](adjust-wps.md) → `adjust_wps/{sample}_wps_out_tss_adjusted.bw` — smoothed and detrended
- [`agg-bw`](agg-bw.md) → `wps/{sample}_wps_out_tss_aggr.wig` — averaged across all TSS windows into one profile

## 4. Controlling it

| What | How |
|------|-----|
| Window extent around each TSS | `--split-interval` |
| Which sites | the `tss` entry of `config/genome.json` for the selected build |
| Sliding window and fragment range | fixed in the workflow at `-W 120`, `120–180 bp` (L-WPS) |
| Threads / memory / walltime | the `wps` entry of `config/cluster.json` |

```json
"wps": {
    "threads": 24,
    "mem": "8G",
    "time": "08:00:00",
    "partition": "norm"
}
```

## 5. Requires

`tss` for the selected genome build. Both bundled builds provide it.

## 6. Reference

Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises an *In Vivo* Nucleosome Footprint that Informs Its Tissues-Of-Origin. *Cell*, 2016;164:57–68.
