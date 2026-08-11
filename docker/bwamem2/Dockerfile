# syntax=docker/dockerfile:1
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Skyler Kuhn <skyler.kuhn@nih.gov>, Ryan Routsong <routsongrm@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="0.0.3"   \
    software="fragmentomics-bwamem2" \
    about.summary="bwa-mem2 + samtools aligner image for the fragmentomics FastQ input path" \
    about.home="https://github.com/OpenOmics/fragmentomics" \
    about.documentation="openomics.github.io/fragmentomics/" \
    about.tags="bwa-mem2 samtools fastqc fastp bedtools alignment fragmentomics"
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone
RUN apt-get update \
    && apt-get -y upgrade
# samtools/bedtools for sort, markdup, proper-pair filtering and indexing;
# fastqc for the read QC that align_fastq runs on both sides of the alignment
# (raw mates, then the analysis BAM) and fastp for read preprocessing; build
# tooling to fetch the pinned bwa-mem2 release binary; python3 to run
# workflow/scripts/python_cpu_arch.py, which picks the bwa-mem2 SIMD variant
# matching the CPU of the node the alignment actually runs on.
RUN apt-get install -y wget git vim samtools bedtools fastqc fastp bzip2 python3-minimal
RUN ln -s /usr/bin/python3 /usr/bin/python
# Pin bwa-mem2 to match the version used to build the reference index
# (config/genome.json -> bwamem2_index). x86_64 prebuilt release binary.
ARG BWAMEM2_VERSION=2.2.1
RUN mkdir /opt2 && cd /opt2 \
    && wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2 \
    && tar -xjf bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2 \
    && mv bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2* /usr/local/bin/ \
    && cd / && rm -rf /opt2
