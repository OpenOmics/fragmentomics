# <code>cfDNAanalyzer <b>feature extraction</b></code>

## 1. About

> **Extracts genomic and fragmentomic feature matrices with a second toolkit, over the same analysis BAMs. Off unless `--cda-features` selects something.**

[cfDNAanalyzer](https://github.com/LiymLab/cfDNAanalyzer) (Zhang *et al.*) is a
separate cfDNA analysis suite from FinaleToolkit, which every other page in this
section documents. It bundles eleven feature extractors — copy number calling
with ichorCNA, nucleosome positioning with DANPOS3, nucleosome profiling with
Griffin, promoter fragmentation entropy with Epic-seq — and this pipeline can
run any subset of them over the same `bams/{sample}.sorted.bam` the
FinaleToolkit analyses read.

Nothing about it runs by default. `--cda-features` names the features to
extract, and an empty selection (the default) leaves cfDNAanalyzer out of the
DAG entirely: no rules are defined and its container is never pulled.

!!! warning "Extraction only. No inference, ever."

    cfDNAanalyzer is three modules: Feature Extraction, Feature
    Processing/Selection, and Machine Learning. **This pipeline runs the first
    one and nothing else.** It passes `--noDA` to the driver, which is precisely
    the cut point: the block that flag guards opens on line 1025 of
    `cfDNAanalyzer.sh` and closes on line 1293, the file's last line, so every
    classifier fit, every cross-validation split, and every per-sample class
    prediction or probability it can produce is inside it.

    No classifier is fit, no feature is selected, no sample is assigned a class
    or a score. The outputs are feature matrices; any inference over them is
    yours to do. `cda_extract` does not take this on trust either — after each
    invocation it asserts that neither `Feature_Processing_and_Selection/` nor
    `Machine_Learning/` was created, and fails the job if one was.

!!! note "The `label` column is a placeholder, not an annotation"

    Every matrix leads with `sample,label`, and `label` is always `0`. This is a
    format artefact: the driver refuses to assemble any matrix without a
    `--labelFile`, and its matrix builder enumerates samples from that file, so
    the workflow synthesizes a one-row label file per sample purely to reach the
    extraction output. Nothing reads it. It is kept in the merged matrices only
    so they stay drop-in comparable with a native cfDNAanalyzer run.

## 2. As the pipeline runs it

Three rules, a scatter and a gather:

```text
cda_regions_bed     once per run, only if a region-specific feature was selected
cda_extract         once per sample
cda_merge_features  once per run, gathers every sample
```

### 2.1 `cda_regions_bed`

Normalizes the regions the region-specific features are scored over into the
sorted BED3 cfDNAanalyzer expects, keeping only the first three columns so any
BED flavor can be handed in:

```text
grep -v -e '^#' -e '^track' -e '^browser' <--cda-regions | tss_interval> \
    | cut -f1,2,3 \
    | sort -k1,1V -k2,2n \
    > cfdnaanalyzer/regions.bed
```

Comment and track lines are dropped rather than passed through, since every
consumer reads this file as plain intervals and would otherwise score a header
line as a region.

### 2.2 `cda_extract`

```text
CDA -I <one-line BAM list> -o <scratch> \
    -F <--cda-features> \
    -g <--genome> \
    -f <reference_fa> \
    -s <pair | single, from the sample's flagstat report> \
    -t <threads> \
    -B <--cda-cna-bin-size> \
    -u <--left-tss-flank> -d <--right-tss-flank> \
    -b cfdnaanalyzer/regions.bed \
    --noDA --labelFile <synthesized>
```

**Why one sample per job.** cfDNAanalyzer walks its `-I` list serially in a
single process, so a single all-samples invocation gets no parallelism at all.
The rule hands it a one-line BAM list instead, and `cda_merge_features` puts the
matrices back together. That is safe because nothing in extraction is
cross-sample: each matrix's columns are derived from an input shared by every
sample (the region BED, the bundled site lists, the 256 end motifs, the fixed
genomic bins), PFE standardization z-scores each row against itself, and the NP
site-list tables are one row per sample.

**Four packaging constraints** shape the rule, and none are optional:

  Sample naming
> cfDNAanalyzer names a sample after its BAM with `.bam` stripped, so the
> analysis BAM is symlinked in as `<sample>.bam`. Passing
> `bams/{sample}.sorted.bam` directly would label every result
> `"{sample}.sorted"`.

---
  The installation tree gets written to, so `CDA` is not invoked in place
> The image puts a `CDA` symlink on `PATH`, and the rule resolves it to find the
> installation rather than hardcoding a path. It cannot just *call* it, though:
> the driver locates itself with `readlink -f "$0"`, which follows the symlink
> all the way back to the real install tree — and that tree is read-only, while
> the `NP` feature generates Griffin snakemake configs and workflows *inside* it.
>
> So the rule builds a shadow install in scratch. The driver is **copied** (a
> symlink would resolve straight back out again), keeping the `CDA` name so logs
> still show what ran; everything beside it is symlinked; and `Griffin/snakemakes`
> is copied and made writable, since that is the part written to.

---
  `NP` runs snakemake itself
> A nested snakemake keeps state in `./.snakemake` and is invoked with
> `--unlock`. The whole thing therefore runs with scratch as its working
> directory — leaving the pipeline's working directory in place would point that
> `--unlock` at the outer workflow's own lock.

---
  `NOF` needs `libR.so` on the loader path
> The rule exports `LD_LIBRARY_PATH="$(R RHOME)/lib:$LD_LIBRARY_PATH"` before
> invoking the driver. `NOF` drives DANPOS3, which reaches R through `rpy2`, and
> `rpy2` embeds R in the Python process rather than going through R's launcher
> script — which is the thing that normally puts `$R_HOME/lib` on the search
> path, by sourcing `$R_HOME/etc/ldpaths`. That directory is on no default path
> otherwise: it has no `/etc/ld.so.conf.d` entry, so `libR.so` is absent from
> the `ldconfig` cache, and the R binary carries no `RUNPATH`.
>
> Embedded, R then starts and fails to `dyn.load` its own base packages, all of
> which link `libR.so`, so `NOF` dies before writing a single wig while the ten
> other features carry on. The image sets this too; the rule sets it as well so
> the feature works against published tags that predate that.

!!! warning "Most features are paired-end only"

    `-s` is not hardcoded — it is read off the `samtools flagstat` report
    `bam_stats` already wrote for the sample, the same way `fastp_bam` does, so a
    single-end BAM is described to cfDNAanalyzer as single-end. That matters
    because seven of the eleven features (`EM`, `FP`, `NP`, `OCF`, `EMR`, `FPR`,
    `PFE`) are paired-end only and cfDNAanalyzer **hard-errors** on them under
    `single`, and because `TSSC` changes how it extends reads based on the flag.

    On single-end data, restrict `--cda-features` to `CNA`, `NOF`, `WPS` and
    `TSSC`. This is a property of cfDNAanalyzer, not of the pipeline, and it is
    not something `--cda-features` validates up front — the layout is only known
    once the BAM has been staged and measured.

!!! warning "cfDNAanalyzer re-filters the BAM, and `--mapscore` does not reach it"

    The driver re-filters every BAM it is given at MAPQ 30, dropping unmapped,
    secondary, QC-fail and duplicate reads (`-q 30 -F 1796`). That is hard-coded
    and not configurable, so these features are measured on a slightly stricter
    subset of the analysis BAM than the FinaleToolkit analyses are.

    Its filter is a floor rather than a ceiling: a `--mapscore` **above** 30
    still holds, because the analysis BAM was already filtered to it before
    cfDNAanalyzer saw it. A `--mapscore` below 30 does not — 30 wins.

### 2.3 `cda_merge_features`

Runs [`workflow/scripts/merge_cda_features.py`](../../workflow/scripts/merge_cda_features.py),
a row concatenation aligned on column name. Where two samples disagree on
columns the union is kept and the mismatch reported, rather than either side
being silently trimmed. A sample that cfDNAanalyzer dropped for failing a
feature's quality control has a header-only CSV and is simply absent from the
merged rows, exactly as it would have been.

## 3. Features

Names are case-insensitive; `all` selects every feature and `none` selects none.

### 3.1 Genome-wide

| Feature | Measures | Matrices |
|---------|----------|----------|
| `CNA` | Copy number alterations, via the bundled ichorCNA, at `--cda-cna-bin-size` | `CNA` |
| `EM` | Fragment end motif frequencies and motif diversity score | `EM_motifs_frequency`, `EM_motifs_mds` |
| `FP` | Short/long fragmentation profile in 100 kb windows | `FP_fragmentation_profile` |

### 3.2 Region-specific

Scored over `--cda-regions`, or the build's `tss_interval` if it is not given.

| Feature | Measures | Matrices |
|---------|----------|----------|
| `NOF` | Nucleosome occupancy and fuzziness (DANPOS3) | `NOF_occupancy`, `NOF_meanfuziness` |
| `NP` | Nucleosome profile at the bundled Griffin TF site lists (Griffin) | `NP_mean_coverage`, `NP_central_coverage`, `NP_amplitude`, plus `NP_site_list/` |
| `WPS` | Windowed protection score, long and short | `WPS_long`, `WPS_short` |
| `OCF` | Orientation-aware cfDNA fragmentation | `OCF` |
| `EMR` | End motif frequencies and MDS, aggregated over and reported per region | `EMR_aggregated_motif_frequency`, `EMR_aggregated_mds`, `EMR_region_motif_frequency`, `EMR_region_mds` |
| `FPR` | Fragmentation profile per region | `FPR_fragmentation_profile_regions` |

### 3.3 Transcription start sites

| Feature | Measures | Matrices |
|---------|----------|----------|
| `PFE` | Promoter fragmentation entropy (Epic-seq) | `PFE` |
| `TSSC` | Average coverage around each TSS, over `--left-tss-flank`/`--right-tss-flank` (deeptools) | `TSSC_average_coverage` |

!!! warning "`WPS` run time scales with the number of regions"

    `WPS` is scored one region at a time by a shell loop. The default region set
    is the build's `tss_interval` — roughly 62,000 regions — which is fine for
    `NOF`, `NP`, `OCF`, `EMR` and `FPR` but makes `WPS` extremely slow. Pair it
    with a focused `--cda-regions` BED.

!!! warning "`PFE` expects a deep targeted panel"

    Epic-seq drops any sample with less than 500x median depth over the
    promoters it scores, and also requires the fragment length mode to fall in
    140–185 bp. A typical whole-genome cfDNA sample fails both, so `PFE.csv` is
    written but may contain only its header. That is not an error and the merge
    tolerates it — the sample is just absent.

## 4. Output

```text
cfdnaanalyzer/
├── regions.bed                     normalized region BED3 (region features only)
├── samples/{sample}/               per-sample matrices, the scatter's output
│   ├── <matrix>.csv
│   └── NP_site_list/               NP only
└── features/                       ← the deliverables
    ├── <matrix>.csv                one row per sample, one column per measurement
    └── NP_site_list/               NP only, one table per TF site list
```

`features/` holds one CSV per matrix in the tables above, each with a leading
`sample,label` pair followed by that feature's measurements. `NP_site_list/`
holds one table per Griffin transcription factor site list — 377 of them for
hg38 — which is why it is tracked as a directory rather than as enumerated
targets: how many there are is a property of the container image.

## 5. Controlling it

| What | How |
|------|-----|
| Which features run (and whether any do) | `--cda-features` |
| Regions the region-specific features score | `--cda-regions` |
| `CNA` bin size | `--cda-cna-bin-size` — 10, 50, 500 or 1000 kb only |
| `TSSC` window | `--left-tss-flank` and `--right-tss-flank`, shared with `cleavage-profile` |
| Threads / memory / walltime | the `cda_extract`, `cda_merge_features` and `cda_regions_bed` entries of `config/cluster.json` |

```json
"cda_extract": {
    "threads": 16,
    "mem": "96G",
    "time": "5-00:00:00",
    "gres": "lscratch:400",
    "partition": "norm"
}
```

The five-day walltime and 400 GB of local scratch are deliberate. `--cda-features all`
runs eleven extractors including a nested Griffin snakemake, and the shadow
install plus intermediate wigs and BEDs all land in `lscratch`. Trim both if you
are selecting only cheap features.

`--cda-cna-bin-size` is restricted because cfDNAanalyzer reads pre-computed GC
and mappability wigs for the requested bin size out of the ichorCNA it bundles,
and those exist only at 10, 50, 500 and 1000 kb. Any other value is rejected by
the frontend rather than failing mid-run.

## 6. Requires

The `cfdnaanalyzer` container image — see
[`docker/cfDNAanalyzer/README.md`](../../docker/cfDNAanalyzer/README.md). Cache
it with `fragmentomics cache` before a run.

cfDNAanalyzer bundles most of its own reference data. Three things come from the
pipeline, and a selected feature that needs a missing one is **silently dropped
from the selection** rather than failing the run:

- `reference_fa` for the selected build — required by `EM`, `EMR` and `NP`
- a region BED — `--cda-regions`, or the build's `tss_interval`, required by
  `NOF`, `NP`, `WPS`, `OCF`, `EMR` and `FPR`
- a build name cfDNAanalyzer accepts, which gates **every** feature rather than
  some, since it is an argument to its driver rather than something one
  extractor reads

Both bundled builds (hg19, hg38) provide all three, so in practice every feature
is available. `hg19` and `hg38` are the only builds cfDNAanalyzer supports at
all: its `CNA`, `FP`, `PFE` and `TSSC` extractors read GC/mappability tracks,
100 kb bin definitions and gene annotations from inside its own installation,
indexed by build name, and its driver rejects any other name outright. A
[`--genome` config file](../pipeline/references.md#5-using-your-own-references)
describing references of your own therefore has to declare which of the two its
coordinates match, with a `cda_genome` entry:

```json
{
    "chrom_sizes":  "/refs/hg19_ucscM.chrom.sizes",
    "reference_fa": "/refs/hg19_ucscM.fa",
    "tss_interval": "/refs/tss.hg19_interval_sorted.bed",
    "cda_genome":   "hg19"
}
```

Without it, `--cda-features` selects nothing and the frontend says so.

## 7. Reference

Zhang J, *et al.* cfDNAanalyzer: a comprehensive toolkit for cell-free DNA
genomic sequencing data analysis. See the
[project repository](https://github.com/LiymLab/cfDNAanalyzer) for the per-feature
method citations (ichorCNA, DANPOS3, Griffin, Epic-seq, DELFI).
