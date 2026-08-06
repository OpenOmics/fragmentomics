# <code>finaletoolkit <b>cleavage-profile</b></code>

## 1. About

> **Calculates cleavage proportion over intervals defined in a BED file.**

Measures, at each base, the proportion of overlapping fragments that *end* there — that is, how often nucleases cut at that exact position. Where [WPS](wps.md) infers nucleosome protection from fragments spanning a window, cleavage profile reads the cut sites directly.

Around a transcription start site the profile shows depleted cleavage where nucleosomes sit and elevated cleavage in the exposed linker regions, giving a complementary view of the same chromatin structure at finer resolution.

## 2. As the pipeline runs it

```text
finaletoolkit cleavage-profile bams/{sample}.sorted.bam <tss> <chrom_sizes> \
    -o cleavage_profile/{sample}_cleavage_profile_tss.bw \
    -l <--left-tss-flank> \
    -r <--right-tss-flank> \
    -q 30 \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -w <threads> \
    -v
```

### 2.1 Positional arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `interval_file`
> **Sites to profile around.**
> *value:* the `tss` file for the selected `--genome` build — individual TSS positions

---
  `chrom_sizes`
> **Chromosome name and length table.**
> *value:* the `chrom_sizes` file for the selected build
>
> Defines the coordinate space of the output bigWig and bounds the flanked intervals so they cannot run off the end of a contig.

### 2.2 Options

  `-l, --left`
> **Base pairs subtracted from each site's start coordinate.**
> *value:* `--left-tss-flank` (default `2000`)
>
> The TSS file gives point positions, so a region has to be built around each one before a profile can be computed. This extends it upstream.

---
  `-r, --right`
> **Base pairs added to each site's stop coordinate.**
> *value:* `--right-tss-flank` (default `2000`)
>
> Extends the region downstream. With both defaults, each profile spans ±2 kb of the TSS — enough to cover the nucleosome-depleted region and the flanking nucleosome array.
>
> !!! note "Asymmetric flanks are allowed"
>
>     `--left-tss-flank` and `--right-tss-flank` are separate options, so the window need not be symmetric. Widening only `--right-tss-flank` extends the profile further into the gene body, where the downstream nucleosome array is most regular.

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `30`

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)
>
> Unlike [`wps`](wps.md), this analysis **does** honor the pipeline's fragment-length options, so the full 50–500 bp range contributes cut sites by default. Passed as `-min`/`-max`; FinaleToolkit also accepts `-lo`/`-hi` as deprecated aliases.

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `24`)

## 3. Output

`cleavage_profile/{sample}_cleavage_profile_tss.bw` — a bigWig of the per-base cleavage proportion across every flanked TSS window.

It is also aggregated across all windows by [`agg-bw`](agg-bw.md) into `cleavage_profile/{sample}_cleavage_profile_aggr.wig`, which is the more interpretable product: averaging over tens of thousands of TSSs is what makes the nucleosome periodicity emerge from what is a sparse signal at any single site.

## 4. Controlling it

| What | How |
|------|-----|
| Window extent | `--left-tss-flank`, `--right-tss-flank` |
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Which sites | the `tss` entry of `config/genome.json` for the selected build |
| Threads / memory / walltime | the `cleavage_profile` entry of `config/cluster.json` |

```json
"cleavage_profile": {
    "threads": 24,
    "mem": "96G",
    "time": "2-00:00:00",
    "partition": "norm"
}
```

This is the most expensive analysis step in the pipeline — 96 GB and a two-day walltime — because it evaluates every fragment endpoint against every base of every flanked TSS window. Widening the flanks increases its cost proportionally.

## 5. Requires

`tss`, `tss_interval` **and** `chrom_sizes` for the selected genome build. (`tss_interval` is needed by the [aggregation](agg-bw.md) step that consumes this output, which is why the workflow gates both on all three.) Both bundled builds provide them.
