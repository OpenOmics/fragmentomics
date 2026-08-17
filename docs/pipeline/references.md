# Reference files

## 1. About

The `--genome` option selects one entry from `config/genome.json`, and that entry supplies every reference file the pipeline uses. Two builds are bundled: `hg38` (the default) and `hg19`.

`--genome` is **required for both input types**. For FastQ input it selects the alignment index in addition to the analysis references; for BAM input it selects the analysis references only, and the alignments are assumed to already match the build. See [Reference contig filter](contig-filter.md) for how that assumption is checked.

All bundled references are built against the **"clean" primary assembly** for their build — the primary chromosomes only, with no alt, decoy or patch scaffolds. Keeping the alignment index, the 2bit sequence, the interval definitions and the chromosome sizes all on the same contig set is what makes fragment coordinates, interval lookups and sequence lookups agree without any per-step coordinate translation.

## 2. What each file is for

  `bwamem2_index`
> **`bwa-mem2` alignment index.**
> *used by:* [`align_fastq`](fastq-alignment.md) — FastQ input only
>
> The pre-built index of the clean primary-assembly FastA. Reads aligned against it can only land on contigs the reference declares, which is why FastQ runs never need the contig filter.
>
> *hg38:* `/data/OpenOmics/references/fragmentomics/hg38/bwamem2_index/hg38_clean.fa`
> *hg19:* `/data/OpenOmics/references/fragmentomics/hg19/bwamem2_index/hg19_clean.fa`

---
  `reference_fa`
> **Clean primary-assembly FastA.**
> *used by:* [`fastq-alignment`](fastq-alignment.md), and the genomic interval BED (below)
>
> The FastA the alignment index, the sequence dictionary and the 2bit file were all derived from.
>
> It is also what the **genomic interval BED is tiled from at run time**. Fixed-width windows tiling the genome used to be built by hand and their path stored in `config/genome.json`; they are now generated per run from this FastA at the width given by [`--interval`](../usage/run.md) (default `1mb`) and written to `intervals/{genome}_{size}_intervals.bed` in the output directory. Those windows are the coordinate space [`frag-length-intervals`](../analyses/frag-length-intervals.md), [`interval-end-motifs`](../analyses/interval-end-motifs.md) and [`delfi`](../analyses/delfi.md) summarize over.
>
> Deriving them rather than bundling them makes the window size a property of the run instead of the build: it is chosen on the command line, recorded in `config.json`, and recoverable from the filename of any result produced with it. It also means this key, not a separate interval key, is what gates the window-based analyses.
>
> *hg38:* `hg38_clean.fa` &nbsp;·&nbsp; *hg19:* `hg19_clean.fa`

---
  `dict`
> **Picard-style sequence dictionary.**
> *used by:* [`filter_reference_contigs`](contig-filter.md) — BAM input only
>
> Supplies the authoritative contig-name and contig-length list for the build. Only the `@SQ` `SN` and `LN` fields are read; `UR` and `M5` are ignored.
>
> *hg38:* `hg38_clean.dict` &nbsp;·&nbsp; *hg19:* `hg19_clean.dict`

---
  `ref2bit`
> **2bit-encoded reference sequence.**
> *used by:* [`end-motifs`](../analyses/end-motifs.md), [`interval-end-motifs`](../analyses/interval-end-motifs.md), [`delfi`](../analyses/delfi.md)
>
> Random-access reference sequence. The motif analyses use it to read the bases at each fragment end; DELFI uses it to compute per-bin GC content for the GC correction.
>
> *hg38:* `hg38.2bit` &nbsp;·&nbsp; *hg19:* `hg19.2bit`

---
  `chrom_sizes`
> **Chromosome name and length table.**
> *used by:* [`delfi`](../analyses/delfi.md), [`adjust-wps`](../analyses/adjust-wps.md), [`cleavage-profile`](../analyses/cleavage-profile.md)
>
> Two-column table defining the coordinate space of the output bigWig files and bounding the intervals DELFI bins over.
>
> *hg38:* `hg38.chrom.sizes` &nbsp;·&nbsp; *hg19:* `hg19.chrom.sizes`

---
---
  `tss`
> **Sorted transcription start site BED.**
> *used by:* [`wps`](../analyses/wps.md), [`cleavage-profile`](../analyses/cleavage-profile.md)
>
> Individual TSS positions. Both analyses build a window around each site and compute a per-base signal across it. The file must be sorted by contig then start, which FinaleToolkit requires.
>
> *hg38:* `tss.hg38_sorted.bed` &nbsp;·&nbsp; *hg19:* `tss.hg19_sorted.bed`

---
  `tss_interval`
> **Sorted TSS interval BED.**
> *used by:* [`coverage`](../analyses/coverage.md), [`adjust-wps`](../analyses/adjust-wps.md), [`agg-bw`](../analyses/agg-bw.md)
>
> The constant-length windows around each TSS, as intervals rather than points. `agg-bw` needs these to know the extent of each window it averages over, and `coverage` uses them as the regions it quantifies.
>
> *hg38:* `tss.hg38_interval_sorted.bed` &nbsp;·&nbsp; *hg19:* `tss.hg19_interval_sorted.bed`

---
  `blacklist`
> **ENCODE blacklist BED.**
> *used by:* [`delfi`](../analyses/delfi.md) (as `--blacklist-file`)
>
> Regions of anomalous mappability or signal. DELFI ignores bins overlapping them, which prevents artefact regions from dominating the short-to-long ratio.
>
> *hg38:* `hg38-blacklist.bed` &nbsp;·&nbsp; *hg19:* `hg19-blacklist.bed`

---
  `gap`
> **UCSC gap track BED4.**
> *used by:* [`delfi`](../analyses/delfi.md) (as `-g`)
>
> Centromere, telomere and short-arm annotations. DELFI excludes these structural gaps from its bins.
>
> *hg38:* `hg38.gap.bed` &nbsp;·&nbsp; *hg19:* `hg19.gap.bed`

## 3. Optional files and graceful degradation

`chrom_sizes`, `reference_fa`, `ref2bit`, `tss` and `tss_interval` gate which analyses run: the workflow only requests an output if the files that step needs are present for the selected build. `gap` and `blacklist` are softer — if either is absent, DELFI runs without the corresponding flag rather than being skipped. `dict` gates the [contig filter](contig-filter.md) alone. See [Pipeline overview §3](overview.md#3-which-steps-run) for the full dependency table.

Both bundled builds define every key, so a default run of either executes the whole workflow.

## 4. Using your own references

On the Biowulf cluster all of the above are already on a shared path and nothing needs to be downloaded. Elsewhere, use <code>fragmentomics <b>install</b></code> to pull the resource bundle, then edit `config/genome.json` so each key points at your local copy.

To add a build, add a new key under `references` with the same file keys. Two constraints are worth stating explicitly:

1. **Every file in a build must describe the same contig set.** Mixing a full-assembly 2bit with a primary-only `reference_fa`, or a `chr`-prefixed TSS file with an unprefixed chrom.sizes, produces empty or misleading output rather than an error. Note that the interval BED inherits its contig set from `reference_fa`, so an assembly FastA carrying alt or decoy scaffolds will tile windows over them too.
2. **Provide a `dict` if you will use BAM input.** Without one, staged BAMs go straight to analysis with no check that they were aligned to your build — a silent failure mode the contig filter exists to prevent.
