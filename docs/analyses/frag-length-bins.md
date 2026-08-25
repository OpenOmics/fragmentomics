# <code>finaletoolkit <b>frag-length-bins</b></code>

## 1. About

> **Retrieves fragment lengths grouped in bins.**

Builds the genome-wide fragment-length distribution for a sample: how many fragments fall into each length bin. This is the most fundamental fragmentomics measurement and the first thing to look at when assessing whether a library is behaving like cfDNA.

Healthy plasma cfDNA shows a characteristic profile — a dominant mononucleosome peak near **167 bp**, a dinucleosome shoulder near **320 bp**, and a ~10 bp periodicity in the sub-nucleosomal range. A sample missing the 167 bp peak, or with a broad distribution extending well past 500 bp, usually indicates contamination with genomic DNA from lysed leukocytes.

This step runs for every genome build; it needs no reference files at all.

## 2. As the pipeline runs it

```text
finaletoolkit frag-length-bins bams/{sample}.sorted.bam \
    -q <--mapscore> \
    --bin-size <--bin-size> \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -p midpoint \
    -o frag_length_bins/{sample}_frag_bin<bin_size>.tsv \
    --histogram-path frag_length_bins/{sample}_frag_bin<bin_size>.png \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `--bin-size`
> **Width of each length bin, in base pairs.**
> *value:* `--bin-size` (default `1`)
>
> At the default of `1`, every distinct fragment length gets its own row — full resolution, which is what preserves the ~10 bp periodicity in the distribution. Larger bins smooth the profile and shrink the output. The value is also embedded in the output filenames, so runs at different bin sizes do not overwrite each other.

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `--mapscore` (default 20)

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)
>
> Fragments outside this window are excluded from the histogram entirely — they are not clamped into the edge bins. Raise `--fragment-maximum` if you need to see the dinucleosome and trinucleosome range in full.

---
  `-p, --intersect-policy`
> **How a fragment is assigned.**
> *value:* `midpoint`
>
> `midpoint` rather than `any` here, so each fragment is counted exactly once. With `any`, a fragment straddling a boundary would contribute to more than one region and distort the counts.

---
  `--histogram-path`
> **Where to write the plotted histogram.**
> *value:* `frag_length_bins/{sample}_frag_bin<bin_size>.png`
>
> A rendered plot of the same distribution as the TSV, produced in the same pass.

## 3. Output

`frag_length_bins/{sample}_frag_bin<bin_size>.tsv` — the binned distribution as a TSV, one row per bin with its fragment count.

`frag_length_bins/{sample}_frag_bin<bin_size>.png` — the same distribution plotted. This is the fastest per-sample QC check in the pipeline; scan it for the 167 bp peak.

!!! note "Summary statistics are not requested"

    FinaleToolkit can append summary statistics as comment lines with `-stats`, and a short-fraction figure with `-sf`. The pipeline does not pass either, so the TSV contains only the bins. Per-interval fragment-length statistics — mean, median, standard deviation, min, max — are available instead from [`frag-length-intervals`](frag-length-intervals.md).

## 4. Controlling it

| What | How |
|------|-----|
| Bin width | `--bin-size` |
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Threads / memory / walltime | the `frag_length_bins` entry of `config/cluster.json` |

```json
"frag_length_bins": {
    "threads": 16,
    "mem": "8G",
    "time": "04:00:00",
    "partition": "norm"
}
```

## 5. Requires

Nothing beyond the analysis BAM. This step runs unconditionally for every sample and every genome build.
