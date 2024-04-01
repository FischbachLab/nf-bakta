FROM mambaorg/micromamba:1.5.7 as app

ARG BAKTA_VER="1.9.3"
ARG DIAMOND_VER="2.1.8"

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="mambaorg/micromamba:1.5.7"
LABEL dockerfile.version="1"
LABEL software="Bakta"
LABEL software.version="${BAKTA_VER}"
LABEL description="rapid & standardized annotation of bacterial genomes, MAGs & plasmids"
LABEL website="https://github.com/oschwengers/bakta"
LABEL license="https://github.com/oschwengers/bakta/blob/main/LICENSE"
LABEL maintainer="Erin Young"
LABEL maintainer.email="eriny@utah.gov"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    procps \
    ca-certificates && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba install --name base -c conda-forge -c bioconda -c defaults bakta=${BAKTA_VER} diamond=${DIAMOND_VER} && \
    micromamba clean -a -y && \
    mkdir /data

ENV PATH="/opt/conda/bin/:$PATH" \
    LC_ALL=C

# Set the working directory
RUN mkdir -p /opt/work
WORKDIR /opt/work

# Set the default command
CMD ["/bin/sh","-c","/bin/bash"]