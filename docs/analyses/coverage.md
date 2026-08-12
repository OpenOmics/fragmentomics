# <code>finaletoolkit <b>coverage</b></code>

## 1. About

> **Calculates fragmentation coverage over intervals defined in a BED file.**

Counts the fragments overlapping each interval and normalizes the counts by the sample's total coverage, so values are comparable across samples with different sequencing depths. In this pipeline the intervals are the TSS windows for the selected genome build, making the output a per-TSS measure of cfDNA abundance.

Coverage is the basis of the project-level [merged workbook](../pipeline/quality-control.md#4-merged-coverage-workbook-merge_coverage_excel), so it is also the primary cross-sample comparison in a run.

## 2. As the pipeline runs it

```text
finaletoolkit coverage bams/{sample}.sorted.bam <tss_interval> \
    -n \
    --scale 1e6 \
    -q <--mapscore> \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -p any \
    -o coverage/{sample}_coverage.bed \
    -w <threads> \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `interval_file`
> **Intervals to quantify coverage over.**
> *value:* the `tss_interval` file for the selected `--genome` build
>
> Constant-length windows centered on transcription start sites. See [Reference files](../pipeline/references.md).

---
  `-n, --normalize`
> **Normalize by total coverage.**
> *set by the pipeline*
>
> Divides each interval's raw fragment count by the sample's genome-wide total, then multiplies by `--scale-factor`. Without this, values would track sequencing depth rather than relative fragment abundance and could not be compared between samples.

---
  `-s, --scale-factor`
> **Scale factor applied to normalized values.**
> *value:* `1e6` (passed as `--scale`)
>
> Turns the normalized fraction into coverage per million fragments, which keeps values in a readable numeric range rather than very small decimals.

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `--mapscore` (default 20)

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)

---
  `-p, --intersect-policy`
> **How a fragment is assigned to an interval.**
> *value:* `any`
>
> A fragment counts toward an interval if it overlaps it anywhere. This is the appropriate choice for coverage — a fragment physically covers part of the interval whether or not its midpoint falls inside. Contrast with [`frag-length-bins`](frag-length-bins.md), which uses `midpoint`.

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `1`)

## 3. Output

`coverage/{sample}_coverage.bed` — a headerless 5-column BED:

```text
chr1	9873	11873	interval_name	14.2731
chr1	28084	30084	interval_name	9.8104
```

| Column | Meaning |
|--------|---------|
| 1 | contig |
| 2 | interval start (0-based) |
| 3 | interval stop |
| 4 | interval name, from the interval BED |
| 5 | normalized, scaled coverage |

These BEDs are merged into `coverage/coverage_summary.xlsx` by the project-level `merge_coverage_excel` step.

## 4. Controlling it

| What | How |
|------|-----|
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Which intervals | the `tss_interval` entry of `config/genome.json` for the selected build |
| Threads / memory / walltime | the `coverage` entry of `config/cluster.json` |

```json
"coverage": {
    "threads": 1,
    "mem": "32000m",
    "time": "72:02:00",
    "gres": "lscratch:200",
    "partition": "norm"
}
```

Coverage is single-threaded by default and is one of the longer-running steps on large BAMs; raise `threads` in `config/cluster.json` to parallelize it.

## 5. Requires

The `tss_interval` file for the selected genome build. Unlike most analyses this step is **not** conditional on that entry being present — it always runs, so a custom build without `tss_interval` fails here rather than skipping the analysis. Both bundled builds provide it.
