# <code>finaletoolkit <b>agg-bw</b></code>

## 1. About

> **Aggregates a bigWig signal over constant-length intervals defined in a BED file.**

Every per-base bigWig the pipeline produces — WPS, adjusted WPS, cleavage profile — covers tens of thousands of separate TSS windows. At any single site the signal is sparse and noisy; the structure only becomes visible when the windows are stacked and averaged.

`agg-bw` does exactly that: it aligns all the intervals on their common coordinate system and averages the signal position by position, producing one profile that represents the sample's average behavior around a TSS. This is the form in which WPS and cleavage results are normally plotted and compared between samples.

## 2. Where the pipeline uses it

The same subcommand runs three times per sample, on three different inputs:

| Workflow step | Input bigWig | Output |
|---------------|-------------|--------|
| `agg_wps` | `wps/{sample}_wps_out_tss.bw` | `wps/{sample}_wps_out_tss_aggr.wig` |
| `agg_adjust_wps` | `adjust_wps/{sample}_wps_out_tss_adjusted.bw` | `adjust_wps/{sample}_wps_out_tss_adj_aggr.wig` |
| `agg_cleavage_profile` | `cleavage_profile/{sample}_cleavage_profile_tss.bw` | `cleavage_profile/{sample}_cleavage_profile_aggr.wig` |

All three are invoked identically apart from input and output paths.

## 3. As the pipeline runs it

```text
finaletoolkit agg-bw <input>.bw <tss_interval> \
    -o <output>.wig \
    --mean \
    -v
```

### 3.1 Arguments

  `input_file`
> **The bigWig to aggregate.**
> *value:* one of the three tracks above

---
  `interval_file`
> **The intervals the signal was computed over.**
> *value:* the `tss_interval` file for the selected `--genome` build
>
> These must be **constant-length** intervals — the aggregation aligns them by position, which is only meaningful if they are all the same size. The bundled `tss_interval` files are the fixed-width windows around each TSS. See [Reference files](../pipeline/references.md).

---
  `-a, --mean`
> **Average instead of taking the median.**
> *set by the pipeline*
>
> Passed as the long form `--mean`. A mean preserves the amplitude of the periodic nucleosome signal linearly, which keeps the aggregate comparable between samples; a median would be more robust to outlier windows but would attenuate the peaks the aggregate exists to reveal. This matches the `--mean` choice in [`adjust-wps`](adjust-wps.md), so the two stages treat the signal consistently.

!!! note "The median filter window is not set"

    `agg-bw` accepts `-m, --median-window-size`, whose help suggests setting it to 120 when aggregating WPS signals. The pipeline does not pass it. For the adjusted-WPS track this is intentional — [`adjust-wps`](adjust-wps.md) has already applied a 200 bp filter, and a second smoothing pass would over-attenuate the nucleosome peaks. The raw-WPS and cleavage aggregates are therefore unsmoothed, so they retain their high-frequency noise; smooth them at plotting time if needed.

## 4. Output

A wiggle (`.wig`) file holding the aggregate signal — one value per position across the common interval extent.

These `_aggr.wig` files are the most immediately usable outputs of the WPS and cleavage analyses. Loaded into a genome browser or plotted directly, the adjusted-WPS aggregate should show the nucleosome-depleted region at the TSS flanked by a regular ~190 bp periodic array.

## 5. Controlling it

| What | How |
|------|-----|
| Which intervals | the `tss_interval` entry of `config/genome.json` for the selected build |
| Aggregation statistic | fixed at `--mean` in the workflow |
| Threads / memory / walltime | the `agg_wps`, `agg_adjust_wps` and `agg_cleavage_profile` entries of `config/cluster.json` |

All three steps share the same modest allocation:

```json
"agg_wps": {
    "threads": 8,
    "mem": "4G",
    "time": "04:00:00",
    "partition": "norm"
}
```

`agg-bw` takes no `-w, --workers` option, so the `threads` value only sizes the job's CPU request.

## 6. Requires

`tss` and `tss_interval` for the aggregation of raw WPS; `chrom_sizes` additionally for the adjusted-WPS and cleavage aggregates, since their upstream steps need it. Both bundled builds provide all three.
