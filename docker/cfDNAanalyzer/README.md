# cfDNAanalyzer image

Feature extraction image for the optional [cfDNAanalyzer](https://github.com/LiymLab/cfDNAanalyzer)
analysis. It carries the whole upstream toolkit, and the `cda_regions_bed`,
`cda_extract` and `cda_merge_features` rules in
[`workflow/rules/cfdnaanalyzer.smk`](../../workflow/rules/cfdnaanalyzer.smk) all
run inside it.

Only cfDNAanalyzer's **Feature Extraction** module is used. The image still
contains its Feature Processing/Selection and Machine Learning modules, because
they ship in the same repository, but the pipeline passes `--noDA` and never
reaches them — see [the analysis page](../../docs/analyses/cfdnaanalyzer.md) for
why that flag is the exact cut point.

Contents:
- **cfDNAanalyzer**, cloned to `/opt/cfDNAanalyzer` and exposed as `CDA` and
  `cfDNAAnalyzer` on `PATH`. `cda_extract` resolves `CDA` to discover the
  installation, so no path into `/opt` is hardcoded in the workflow and moving
  the install only means moving the symlink. It does not *invoke* `CDA` in
  place: the driver locates itself with `readlink -f "$0"`, which follows the
  symlink back to this read-only tree, and the `NP` feature writes generated
  Griffin workflows into `Griffin/snakemakes` underneath it. The driver is
  therefore copied into scratch as `CDA` and run from a shadow install
- **R 4.3.x / Bioconductor 3.18**, from the `bioconductor/bioconductor_docker`
  base image. The release is not arbitrary — it matches the
  `BiocManager::install(version = "3.18")` pin in upstream's
  `install_R_packages.R`. The base image also supplies the C/C++/Fortran
  toolchain and system `-dev` libraries those packages build against, plus a
  binary package repo so most install prebuilt
- **Python 3.7.16** and the bioinformatics stack, in a micromamba env created
  from upstream's `environment.yml` (`samtools`, `bedtools`, `pysam`,
  `snakemake` 5.19.2, `pandas`, `scikit-learn`, ...). The env is first on `PATH`;
  R deliberately is not in it and stays on the base image's `PATH`
- **`jq`**, from apt. Not optional and not a build tool: `cfDNAanalyzer.sh`
  checks for it before doing anything else and exits if it is missing, so
  without it *every* invocation fails whichever features were requested. It is
  what subsets `Feature_Processing/config_all_feature.json` down to the selected
  features, which is how the matrix builder learns what to assemble
- **`deeptools` 3.5.3**, pip-installed into the env. Upstream's
  `environment.yml` pins `deeptoolsintervals` (a *dependency* of deeptools) but
  not deeptools itself, so `bamCoverage` and `multiBigwigSummary` are otherwise
  absent and the `TSSC` feature fails at its first step. 3.5.3 is the newest
  release still declaring Python 3.7 support; it is installed `--no-deps`
  because all ten of its requirements are already satisfied by
  `environment.yml`'s exact-build pins and letting pip resolve them would churn
  the pinned env
- **`rpy2` 3.3.3 in ABI mode.** `environment.yml` pip-builds rpy2 during env
  creation, and there is no R inside the conda env for it to compile against.
  `RPY2_CFFI_MODE=ABI` at build time makes it bind to `libR` at run time
  through cffi instead, and `R_HOME=/usr/local/lib/R` is what points it at the
  base image's R

## Executable bits

Upstream commits **every** file in the repository mode `100644` — not one
carries the executable bit — but `cfDNAanalyzer.sh` runs several of the bundled
helpers as commands rather than through an interpreter. The Dockerfile marks
each of those explicitly, asserting the path exists first so an upstream move or
rename fails the build instead of failing later inside a feature:

| Feature | Helper | Why it is executed directly |
|---------|--------|-----------------------------|
| `CNA`  | `ichorCNA/hmmcopy_utils/bin/readCounter` | bins the BAM into a wig |
| `NOF`  | `DANPOS3/wigToBigWig` | nucleosome wig to bigWig |
| `NOF`  | `DANPOS3/bigWigAverageOverBed` | occupancy per region |
| `WPS`  | `WPS/calculate_wps.sh` | per-region driver |
| `WPS`  | `WPS/samtools` | bundled build, needs `-m`/`-M` fragment length filters |
| `WPS`  | `WPS/FilterUniqueBAM.py` | piped, resolved by shebang |
| `WPS`  | `WPS/extractReadStartsFromBAM2Wig.py` | piped, resolved by shebang |
| `OCF`  | `OCF/OCF.sh` | per-region driver |
| `OCF`  | `OCF/bedtools` | bundled build, called by `OCF.sh` |
| `PFE`  | `Epic-seq/code/epic_wrapper.sh` | invoked by `runEPIC.R` |

Note that the `chmod -R a+rX /opt/cfDNAanalyzer` in the same stage does **not**
cover these. Capital `X` only adds `+x` to directories and to files that are
already executable, which in this tree is none of them.

## Reference files

cfDNAanalyzer bundles its own references, so nothing under
`/data/OpenOmics/references/fragmentomics/` is needed for most features and
`config/genome.json` gains no new keys. What ships in the image:

```
ichorCNA/ichorCNA/inst/extdata/   GC and mappability wigs, hg19 + hg38, at
                                  10/50/500/1000 kb — which is why
                                  --cda-cna-bin-size accepts only those four
Fragmentation_profile/            {hg19,hg38}_100kb_WG.txt genome-wide windows
DANPOS3/                          {hg19,hg38}.chrom.sizes
Griffin/Ref_{hg19,hg38}/sites/    transcription factor site lists for NP
TSScoverage/                      TSS_{hg19,hg38}_uniq.bed
Epic-seq/code/priordata/          PFE control and TSS prior data
```

Two inputs do come from the pipeline rather than the image: the reference FastA
(`reference_fa`, required by `EM`, `EMR` and `NP`) and the region BED that the
region-specific features are scored over, which is `--cda-regions` when given
and the build's `tss_interval` otherwise.

## Build & push

```bash
docker build --platform linux/amd64 -t rroutsong/fragmentomics_cfdnaanalyzer:0.0.1 .
docker push rroutsong/fragmentomics_cfdnaanalyzer:0.0.1
```

The image URI is registered in `config/containers.json` under the
`cfdnaanalyzer` key and can be cached locally with `fragmentomics cache`.

This is a large image — a Bioconductor base, a full R package set, a conda env
and the bundled reference data — so expect a long first build and a slow initial
pull. It is only pulled when a run actually selects `--cda-features`.
