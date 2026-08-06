# <code>finaletoolkit <b>delfi</b></code>

## 1. About

> **Calculates DELFI features over the genome, returning information about (GC-corrected) short fragments, long fragments, DELFI ratio, and total fragments.**

DELFI — **DE**tection of **F**ragmentation for **L**iquid biopsy **I**nvestigation, from Cristiano *et al.* (2019) — measures the ratio of short to long cfDNA fragments in bins across the genome. Tumor-derived cfDNA is systematically shorter than cfDNA from healthy cells, so the *spatial pattern* of this ratio is one of the strongest published cfDNA cancer signals.

The measurement is GC-corrected, which is essential: bin-level GC content confounds both fragment length and coverage, so an uncorrected ratio partly reports GC bias rather than biology. DELFI is the most reference-hungry analysis in the pipeline, taking four positional files plus the gap and blacklist tracks.

## 2. As the pipeline runs it

```text
finaletoolkit delfi bams/{sample}.sorted.bam <chrom_sizes> <ref2bit> <intervals> \
    -q 30 \
    --blacklist-file <blacklist> \
    -g <gap> \
    -o delfi/{sample}_delfi.bed \
    -w <threads> \
    -v \
    --no-merge-bins
```

### 2.1 Positional arguments

The four positionals must be given in this order:

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `chrom_sizes`
> **Chromosome name and length table.**
> *value:* the `chrom_sizes` file for the selected `--genome` build
>
> Bounds the region DELFI bins over.
>
> !!! note
>
>     FinaleToolkit's help notes that to replicate the original Cristiano *et al.* methodology this file "should contain only autosomes". The bundled `chrom_sizes` files are full builds including sex chromosomes and chrM, so DELFI bins here span more than the original publication's autosome-only scope. Filter the file if you need strict replication.

---
  `reference_file`
> **2bit reference sequence.**
> *value:* the `ref2bit` file for the selected build
>
> Used to compute per-bin GC content for the GC correction.

---
  `bins_file`
> **Bins to compute DELFI over.**
> *value:* the `intervals` file for the selected build — 5 kb windows
>
> !!! warning "The bundled bins are 5 kb, not 100 kb"
>
>     The original DELFI methodology uses 100 kb bins merged to 5 Mb. This pipeline passes its general-purpose 5 kb interval file, giving finer bins than the publication. Combined with `--no-merge-bins` below, the output is at 5 kb resolution throughout. This is a deliberate choice — it keeps a single interval definition across `frag-length-intervals`, `interval-end-motifs` and `delfi` — but it means the values are not numerically comparable to published DELFI scores without re-binning.

### 2.2 Options

  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `30`

---
  `-b, --blacklist-file`
> **Regions to ignore.**
> *value:* the `blacklist` entry for the selected build, passed as `--blacklist-file`
>
> ENCODE blacklist regions — anomalous mappability or signal artefacts. Excluding them prevents a handful of pathological loci from dominating the ratio. If the selected build has no `blacklist` entry, the flag is omitted and DELFI runs without it.

---
  `-g, --gap-file`
> **Structural gaps to exclude.**
> *value:* the `gap` entry for the selected build
>
> A BED4 of centromere, telomere and short-arm annotations from the UCSC gap track. Also omitted if the build does not define it.

---
  `-M, --no-merge-bins`
> **Keep the input bins; do not merge to 5 Mb.**
> *set by the pipeline*
>
> By default DELFI merges its bins up to `--window-size` (5 Mb) before reporting. The pipeline suppresses that so the output stays at the resolution of the input interval file, leaving any aggregation to downstream analysis rather than baking it in.

---
  `-w, --workers`
> **Worker processes.**
> *value:* the step's `threads` from `config/cluster.json` (default `24`)

!!! note "GC correction is on; the hg19 no-coverage fix is not disabled"

    The pipeline passes neither `-G, --no-gc-correct` (so GC correction **is** applied) nor `-R, --keep-nocov`. `-R` exists to skip the removal of two hg19-specific no-coverage regions and FinaleToolkit's help advises setting it for non-hg19 references; the pipeline leaves it unset for both builds, so those two regions are removed regardless of build. The affected span is negligible relative to a genome-wide bin set.

## 3. Output

`delfi/{sample}_delfi.bed` — one row per bin with its short-fragment count, long-fragment count, DELFI ratio and total fragment count, GC-corrected.

## 4. Controlling it

| What | How |
|------|-----|
| Bin resolution | the `intervals` entry of `config/genome.json` for the selected build |
| Blacklist / gap exclusion | the `blacklist` and `gap` entries — remove a key to drop the flag |
| Threads / memory / walltime | the `delfi` entry of `config/cluster.json` |

```json
"delfi": {
    "threads": 24,
    "mem": "16G",
    "time": "08:00:00",
    "partition": "norm"
}
```

Note that DELFI does not take `--fragment-minimum`/`--fragment-maximum` — it defines short and long internally, since that partition *is* the measurement.

## 5. Requires

`ref2bit`, `intervals` **and** `chrom_sizes` for the selected genome build. `blacklist` and `gap` are optional; each is passed only if present. Both bundled builds provide all five.

## 6. Reference

Cristiano S, *et al.* Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature*, 2019;570:385–389.
