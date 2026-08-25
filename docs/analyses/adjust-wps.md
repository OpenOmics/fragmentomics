# <code>finaletoolkit <b>adjust-wps</b></code>

## 1. About

> **Adjusts raw Windowed Protection Score (WPS) by applying a median filter and Savitsky-Golay filter.**

Raw [WPS](wps.md) is not directly interpretable. It carries two artefacts that swamp the nucleosome signal:

1. **A long-wavelength trend** driven by local coverage — regions with more fragments have systematically higher scores regardless of nucleosome occupancy.
2. **High-frequency noise** from the finite number of fragments contributing to any single base.

This step removes both. A moving-average filter over a 200 bp window subtracts the local baseline, and edge subtraction removes the residual offset at each window's flanks. The result is a detrended track in which the ~190 bp nucleosome periodicity around a TSS is visible.

This step reads a bigWig, not the BAM, so it does not touch fragment data.

## 2. As the pipeline runs it

```text
finaletoolkit adjust-wps wps/{sample}_wps_out_tss.bw <tss_interval> <chrom_sizes> \
    -o adjust_wps/{sample}_wps_out_tss_adjusted.bw \
    -i <--split-interval> \
    -m 200 \
    -S \
    --subtract-edges \
    -v
```

### 2.1 Positional arguments

  `input_file`
> **Raw WPS bigWig.**
> *value:* `wps/{sample}_wps_out_tss.bw` from [`wps`](wps.md)

---
  `interval_file`
> **The intervals WPS was computed over.**
> *value:* the `tss_interval` file for the selected `--genome` build
>
> The TSS windows as *intervals*, so the filter knows the extent of each region it is detrending. Note that [`wps`](wps.md) upstream takes the `tss` point file, while this step takes the `tss_interval` file — see [Reference files](../pipeline/references.md).

---
  `chrom_sizes`
> **Chromosome name and length table.**
> *value:* the `chrom_sizes` file for the selected build
>
> Defines the coordinate space of the output bigWig.

### 2.2 Options

  `-i, --interval_size`
> **Size of each interval in the interval file.**
> *value:* `--split-interval` (default `5000`)
>
> Must match the `-i` given to [`wps`](wps.md), which the pipeline guarantees by passing the same option to both.

---
  `-m, --median-window-size`
> **Filter window, in base pairs.**
> *value:* `200`
>
> The window whose average is subtracted from each position. At 200 bp it is slightly wider than a nucleosome footprint, so it captures the local baseline without flattening the nucleosome peaks themselves.

---
  `-S, --exclude-savgol`
> **Skip Savitzky–Golay filtering.**
> *set by the pipeline*
>
> The polynomial smoothing pass is disabled, leaving the moving-average detrending as the only filter. Savitzky–Golay smoothing would attenuate the sharp nucleosome peaks the analysis is meant to resolve. With `-S` set, the `-s, --savgol-window-size` and `-p, --savgol-poly-deg` options have no effect.

---
  `--mean`
> **Use a mean filter instead of a median filter.**
> *set by the pipeline*
>
> Despite the option being named `--median-window-size`, `--mean` makes the filter a moving average. A mean filter responds linearly to the underlying signal, which is the appropriate behavior for subtracting a coverage-driven baseline; a median filter is more robust to outliers but distorts the amplitude of the periodic signal being preserved.

---
  `--subtract-edges`
> **Subtract the median of each window's flanks from the whole interval.**
> *set by the pipeline*
>
> Takes the median of the first and last 500 bases of each interval and subtracts it, which removes the per-window offset remaining after detrending. This is what makes adjusted values comparable *between* TSS windows, and therefore what makes the [aggregate](agg-bw.md) meaningful. The 500 bp flank size is FinaleToolkit's default (`--edge-size`); the pipeline does not override it.

!!! note "Single-process step"

    `adjust-wps` accepts `-w, --workers`, but the pipeline does not pass it — the work is a filter pass over a bigWig rather than a scan over reads. Its `threads` allocation in `config/cluster.json` still applies to the job's CPU request.

## 3. Output

`adjust_wps/{sample}_wps_out_tss_adjusted.bw` — the detrended, edge-corrected per-base WPS track.

This is the WPS product to interpret. It is also aggregated across all TSS windows by [`agg-bw`](agg-bw.md) into `adjust_wps/{sample}_wps_out_tss_adj_aggr.wig`.

## 4. Controlling it

| What | How |
|------|-----|
| Interval size | `--split-interval` (must match [`wps`](wps.md)) |
| Which intervals | the `tss_interval` entry of `config/genome.json` for the selected build |
| Filter window, filter type, edge subtraction | fixed in the workflow at `-m 200`, `--mean`, `--subtract-edges` |
| Threads / memory / walltime | the `adjust_wps` entry of `config/cluster.json` |

```json
"adjust_wps": {
    "threads": 24,
    "mem": "32G",
    "time": "1-00:00:00",
    "partition": "norm"
}
```

The 24-hour walltime reflects that this is a per-base pass over every TSS window in the genome.

## 5. Requires

`tss`, `tss_interval` **and** `chrom_sizes` for the selected genome build, plus the raw WPS bigWig from [`wps`](wps.md). Both bundled builds provide all three files.
