# =============================================================================
# cfDNAanalyzer  (LiymLab/cfDNAanalyzer)
# https://github.com/LiymLab/cfDNAanalyzer
#
# Base: official Bioconductor image (Ubuntu). It ships R 4.3.x, Bioconductor
# 3.18 (matching the `BiocManager::install(version = "3.18")` pin in
# install_R_packages.R), a full C/C++/Fortran build toolchain, every system
# -dev library the R packages need, and a binary package repo so most R
# packages install prebuilt instead of compiling from source.
#
# The Python/bioinformatics stack (Python 3.7.16, samtools, bedtools, ...) is
# provided by a micromamba env created from environment.yml.
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_18

ARG DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# -----------------------------------------------------------------------------
# Only the bits micromamba/clone need; R build deps already ship in the image.
#
# jq is not optional and is not a build tool: cfDNAanalyzer.sh checks for it
# before doing anything else and exits if it is absent, so without it *every*
# invocation fails regardless of which features were requested. It uses jq to
# subset Feature_Processing/config_all_feature.json down to the selected
# features, which is what tells the matrix builder what to assemble.
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates curl bzip2 xz-utils unzip procps jq \
 && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Install micromamba (single static binary, ~10 MB). No Miniconda, no user.
# -----------------------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

RUN mkdir -p /opt/conda \
 && curl -fSL --retry 5 --retry-delay 3 --retry-all-errors \
        -o /tmp/micromamba.tar.bz2 \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
 && tar -xvj -C /usr/local -f /tmp/micromamba.tar.bz2 bin/micromamba \
 && rm -f /tmp/micromamba.tar.bz2 \
 && micromamba shell init -s bash -r /opt/conda >/dev/null

# -----------------------------------------------------------------------------
# Clone the repo.
# -----------------------------------------------------------------------------
WORKDIR /opt
RUN git clone --depth 1 https://github.com/LiymLab/cfDNAanalyzer.git \
 && chmod -R a+rX /opt/cfDNAanalyzer

WORKDIR /opt/cfDNAanalyzer

# -----------------------------------------------------------------------------
# Python/bioinformatics env from environment.yml (name: cfDNAanalyzer):
# Python 3.7.16, samtools, bedtools, deeptools, etc. environment.yml does NOT
# install R; R comes from the base image.
#
# RPY2_CFFI_MODE=ABI: environment.yml pip-builds rpy2==3.3.3 during env create
# and there is no R inside this conda env, so ABI mode makes rpy2 bind to libR
# at runtime (cffi) instead of requiring R headers at build time. R_HOME (set
# below) points it at the image's R.
# -----------------------------------------------------------------------------
RUN RPY2_CFFI_MODE=ABI micromamba env create -y -f /opt/cfDNAanalyzer/environment.yml \
 && micromamba clean --all --yes

# -----------------------------------------------------------------------------
# deeptools, which environment.yml omits: it pins `deeptoolsintervals` (a
# dependency of deeptools) but not deeptools itself, so `bamCoverage` and
# `multiBigwigSummary` are absent and the TSSC feature fails at its first step.
#
# 3.5.3 is the newest release that still declares support for Python 3.7, which
# is what environment.yml pins the env to. --no-deps is deliberate: all ten of
# its requirements (numpy, scipy, matplotlib, pysam, numpydoc, pyBigWig, py2bit,
# plotly, deeptoolsintervals, importlib-metadata) are already satisfied by
# environment.yml's exact-build pins, and letting pip resolve them would churn
# the pinned env.
# -----------------------------------------------------------------------------
RUN /opt/conda/envs/cfDNAanalyzer/bin/pip install --no-cache-dir --no-deps \
        deeptools==3.5.3 \
 && /opt/conda/envs/cfDNAanalyzer/bin/bamCoverage --version \
 && /opt/conda/envs/cfDNAanalyzer/bin/multiBigwigSummary --version

# -----------------------------------------------------------------------------
# PATH: conda env bin first (python, samtools, ...); the image's R stays on the
# default PATH (/usr/local/bin) since the conda env has no R of its own.
# R_HOME lets rpy2 (ABI mode) locate libR at runtime.
#
# LD_LIBRARY_PATH is needed too, and only by rpy2. $R_HOME/lib is on no default
# loader path -- no /etc/ld.so.conf.d entry, no RUNPATH in the binary -- and
# what normally covers that is R's launcher script, which sources
# $R_HOME/etc/ldpaths. rpy2 embeds R rather than running the launcher, so R
# starts and then cannot dyn.load its own base packages, all of which link
# libR.so. DANPOS3 is the only rpy2 consumer here and the NOF feature its only
# caller, so without this NOF is the one feature that cannot run. That directory
# holds nothing but libR.so, so putting it on the search path shadows nothing.
# -----------------------------------------------------------------------------
ENV CONDA_PREFIX=/opt/conda/envs/cfDNAanalyzer
ENV R_HOME=/usr/local/lib/R
ENV LD_LIBRARY_PATH=/usr/local/lib/R/lib:$LD_LIBRARY_PATH
ENV PATH=/opt/conda/envs/cfDNAanalyzer/bin:$PATH

# -----------------------------------------------------------------------------
# Install the CRAN + Bioconductor R packages, using the image's R 4.3 / Bioc
# 3.18. Runs with an explicitly clean PATH-resolved R (base image R), not the
# conda env.
# -----------------------------------------------------------------------------
RUN Rscript /opt/cfDNAanalyzer/install_R_packages.R

# -----------------------------------------------------------------------------
# Assert rpy2 can embed R, the operation the LD_LIBRARY_PATH above exists for.
# It regresses easily and presents as an NOF bug hours into a run rather than as
# an image bug, so it is caught here instead. danpos.py runs from its own
# directory (its modules import each other by bare name) and with no arguments
# prints a version banner, which is enough to prove the import chain resolves.
# -----------------------------------------------------------------------------
RUN /opt/conda/envs/cfDNAanalyzer/bin/python -c \
      "import rpy2.robjects as ro; print('rpy2 embeds', ro.r('R.version.string')[0])" \
 && cd /opt/cfDNAanalyzer/DANPOS3 \
 && /opt/conda/envs/cfDNAanalyzer/bin/python danpos.py | grep -q 'danpos version'

# -----------------------------------------------------------------------------
# cfDNAAnalyzer execution setup
# -----------------------------------------------------------------------------
RUN chmod -R a+rX /opt/conda /opt/cfDNAanalyzer \
 && chmod +x /opt/cfDNAanalyzer/cfDNAanalyzer.sh \
 && ln -s /opt/cfDNAanalyzer/cfDNAanalyzer.sh /usr/bin/cfDNAAnalyzer \
 && ln -s /opt/cfDNAanalyzer/cfDNAanalyzer.sh /usr/bin/CDA

# -----------------------------------------------------------------------------
# Executable bits for the bundled helpers cfDNAanalyzer.sh runs as commands
# rather than through an interpreter.
#
# Upstream commits every file in the repo mode 100644 -- not one carries the
# executable bit -- so each of these has to be marked here or the feature that
# calls it dies with "Permission denied". The `chmod -R a+rX` above does not
# cover them: capital X only adds +x to directories and to files that are
# already executable, which in this tree is none of them.
#
#   CNA   ichorCNA/hmmcopy_utils/bin/readCounter   bins the BAM into a wig
#   NOF   DANPOS3/wigToBigWig                      nucleosome wig to bigWig
#   NOF   DANPOS3/bigWigAverageOverBed             occupancy per region
#   WPS   WPS/calculate_wps.sh                     per-region driver
#   WPS   WPS/samtools                             bundled; needs -m/-M support
#   WPS   WPS/FilterUniqueBAM.py                   piped, resolved by shebang
#   WPS   WPS/extractReadStartsFromBAM2Wig.py      piped, resolved by shebang
#   OCF   OCF/OCF.sh                               per-region driver
#   OCF   OCF/bedtools                             bundled; called by OCF.sh
#   PFE   Epic-seq/code/epic_wrapper.sh            invoked by runEPIC.R
#
# Each path is asserted before being chmod'ed, so an upstream change that moves
# or renames one of them fails this build rather than failing at runtime inside
# whichever feature needed it.
# -----------------------------------------------------------------------------
RUN set -eux; \
    for tool in \
        ichorCNA/hmmcopy_utils/bin/readCounter \
        DANPOS3/wigToBigWig \
        DANPOS3/bigWigAverageOverBed \
        WPS/calculate_wps.sh \
        WPS/samtools \
        WPS/FilterUniqueBAM.py \
        WPS/extractReadStartsFromBAM2Wig.py \
        OCF/OCF.sh \
        OCF/bedtools \
        Epic-seq/code/epic_wrapper.sh \
    ; do \
        test -f "/opt/cfDNAanalyzer/${tool}"; \
        chmod a+x "/opt/cfDNAanalyzer/${tool}"; \
    done

# -----------------------------------------------------------------------------
# provenance
# -----------------------------------------------------------------------------
ADD Dockerfile /Dockerfile
RUN { \
      echo "=== SYSTEM ==="; uname -a; cat /etc/os-release; \
      echo; echo "=== ENV ==="; env | sort; \
      echo; echo "=== APT PACKAGES ==="; dpkg-query -W -f='${Package}\t${Version}\n' 2>/dev/null || echo "dpkg not available"; \
      echo; echo "=== PYTHON ==="; python3 --version 2>&1; pip3 list --format=freeze 2>/dev/null || echo "pip not available"; \
      echo; echo "=== R PACKAGES ==="; R --version 2>&1 | head -1; Rscript -e 'ip <- installed.packages()[,c("Package","Version")]; write.table(ip, quote=FALSE, row.names=FALSE, sep="\t")' 2>/dev/null || echo "R not available"; \
    } > /session.txt 2>&1

# -----------------------------------------------------------------------------
# No ENTRYPOINT (Singularity ignores/overrides ENTRYPOINT behavior anyway).
# Default CMD is a shell; users invoke cfDNAanalyzer.sh directly.
# -----------------------------------------------------------------------------
WORKDIR /work
CMD ["bash"]
