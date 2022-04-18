FROM ubuntu:22.04

#Metadata
LABEL base.image="ubuntu:22.04"
LABEL version="0.2"
LABEL description="TimeAttackGenComp tools, including SNAP 2.0.1 aligner and BCFtools 1.15.1"
LABEL tag="TimeAttackGenComp"

RUN buildDeps='wget gcc g++ make bzip2' && \
    apt-get update -qq && \
    apt-get upgrade -y && \
    apt-get clean all && \
    apt-get install -y --no-install-recommends \
        r-base r-cran-rcolorbrewer liblzma-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev && \
    apt-get install -y --no-install-recommends \
        $buildDeps && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* && \
    wget https://github.com/teerjk/TimeAttackGenComp/archive/refs/heads/master.zip && \
    unzip master.zip && \
    wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2 && \
    tar -xvjf bcftools-1.15.1.tar.bz2 && \
    cd bcftools-1.15.1/ && \
    make && cd ../ && \
    wget https://github.com/amplab/snap/archive/refs/tags/v2.0.1.tar.gz && \
    tar -xvzf v2.0.1.tar.gz && \
    cd snap-2.0.1/ && \
    make && cd ../ && \
    rm master.zip bcftools-1.15.1.tar.bz2 v2.0.1.tar.gz && \
    apt-get purge -y --auto-remove $buildDeps
