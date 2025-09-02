# syntax=docker/dockerfile:1
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Skyler Kuhn <skyler.kuhn@nih.gov>, Ryan Routsong <routsongrm@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="0.0.1"   \
    software="fragmentomics" \
    about.summary="An awesome fragmentomics pipeline" \
    about.home="https://github.com/OpenOmics/fragmentomics" \
    about.documentation="openomics.github.io/fragmentomics/" \
    about.tags="finaltools suite for fragmentomics"
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone
RUN apt-get update \
    && apt-get -y upgrade
RUN apt-get install -y wget git vim samtools bedtools python3-pip golang-go
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN mkdir /opt2; cd /opt2; wget https://github.com/shenwei356/seqkit/releases/download/v2.10.1/seqkit_linux_arm64.tar.gz
RUN cd /opt2; tar -xzvf seqkit_linux_arm64.tar.gz; mv seqkit /usr/bin; cd /; rm -rf /opt2
RUN pip install finaletoolkit
