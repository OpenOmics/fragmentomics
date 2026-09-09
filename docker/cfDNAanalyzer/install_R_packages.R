#!/usr/bin/env Rscript
#
# R packages for cfDNAanalyzer. A rewrite of upstream's install_R_packages.R
# (https://github.com/LiymLab/cfDNAanalyzer/blob/main/install_R_packages.R),
# with the same package sets and the same pinned versions. Three differences:
#
#   1. Mirror. Upstream hardcodes Tsinghua University's, the right pick for its
#      authors and a poor one from us-east: 2.7 s to first byte on an archived
#      tarball, against 0.12 s for the first candidate below. Candidates are
#      tried in order and the first that serves an index wins, so a mirror that
#      is down or blocked on the build network does not fail the build.
#   2. Only the CRAN entry is set, where upstream replaces `repos` wholesale.
#      Bioconductor therefore keeps resolving through BioCcontainers, the
#      prebuilt binary repo the base image configures, and those packages
#      install without compiling.
#   3. install.packages reports a failed install as a warning rather than an
#      error, so R exits 0 having left packages missing and the feature that
#      needed one dies hours into a run instead. The check at the end turns that
#      back into a nonzero exit.
#
# The pinned CRAN versions install from source tarballs in the CRAN archive,
# which is the only place install_version can satisfy an exact version, so those
# build from source whatever the mirror.

MIRRORS <- c(
    "https://cloud.r-project.org",                     # CRAN's CDN, Ashburn VA edge
    "https://archive.linux.duke.edu/cran",             # Durham NC
    "https://p3m.dev/cran/__linux__/jammy/2024-04-23"  # Posit; base image default
)

CRAN_PINNED <- c(
    DescTools   = "0.99.40",
    zoo         = "1.8.12",
    plyr        = "1.8.9",
    reshape2    = "1.4.4",
    data.table  = "1.15.2",
    MASS        = "7.3-60.0.1",
    e1071       = "1.7-14",
    gtools      = "3.9.5",
    matrixStats = "1.2.0",
    optparse    = "1.7.4",
    httr        = "1.4.7",
    tidyverse   = "2.0.0",
    RCurl       = "1.98-1.14",
    devtools    = "2.4.5"
)

BIOC_VERSION <- "3.18"

BIOCONDUCTOR <- c(
    "HMMcopy",
    "GenomeInfoDb",
    "GenomicRanges",
    "Rsamtools",
    "GenomicAlignments",
    "biovizBase",
    "rtracklayer",
    "BSgenome.Hsapiens.UCSC.hg19",
    "BSgenome.Hsapiens.UCSC.hg38"
)

# Named by the package the repo installs, so the closing check can look for it.
GITHUB <- c(ichorCNA = "broadinstitute/ichorCNA")


serves_an_index <- function(mirror) {
    isTRUE(tryCatch(
        nrow(suppressWarnings(
            available.packages(contrib.url(mirror, "source"))
        )) > 0,
        error = function(e) FALSE
    ))
}

first_reachable <- function(mirrors) {
    for (mirror in mirrors) {
        if (serves_an_index(mirror)) {
            return(mirror)
        }
        message("CRAN mirror unreachable, trying next: ", mirror)
    }
    stop("no CRAN mirror reachable; check egress from the build network")
}


cran <- first_reachable(MIRRORS)
message("CRAN mirror: ", cran)
options(repos = c(CRAN = cran))

for (pkg in c("remotes", "BiocManager")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

for (pkg in names(CRAN_PINNED)) {
    remotes::install_version(pkg, version = CRAN_PINNED[[pkg]])
}

BiocManager::install(version = BIOC_VERSION)
BiocManager::install(BIOCONDUCTOR)

for (repo in GITHUB) {
    remotes::install_github(repo)
}

wanted <- c(names(CRAN_PINNED), BIOCONDUCTOR, names(GITHUB))
missing <- setdiff(wanted, rownames(installed.packages()))
if (length(missing)) {
    stop("R packages missing after install: ", paste(missing, collapse = ", "))
}
message("all ", length(wanted), " requested R packages present")
