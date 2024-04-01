# Setup base image
FROM mambaorg/micromamba:1.5.8 as app

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    procps \
    ca-certificates && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*
RUN  micromamba install --name base -c conda-forge -c bioconda -c default -y python=3.10 numpy polars pydantic pyhmmer \
    && micromamba clean --all --yes \
    && mkdir -p /opt/work

ENV PATH="/opt/conda/bin/:$PATH" \
    LC_ALL=C


WORKDIR /opt/work

# Metadata
LABEL container.maintainer="Sunit Jain" \
    container.base.image="continuumio/miniconda3" \
    software.name="pyhmmer" \
    software.description="Cython bindings and Python interface to HMMER3." \
    software.website="https://github.com/althonos/pyhmmer" \
    container.category="aligner"

RUN chmod -R +rx ./
# Set the default command
CMD ["/bin/sh","-c","/bin/bash"]