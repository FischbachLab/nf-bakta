ARG HMMER_VER="3.4"

## Builder ##
FROM ubuntu:jammy as builder

ARG HMMER_VER

RUN apt-get update && apt-get install --no-install-recommends -y \
    wget \
    build-essential &&\
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN wget http://eddylab.org/software/hmmer/hmmer-${HMMER_VER}.tar.gz &&\
    tar -xvf hmmer-${HMMER_VER}.tar.gz &&\
    cd hmmer-${HMMER_VER} &&\
    ./configure &&\
    make && make install

## App ##
FROM ubuntu:jammy as app

ARG HMMER_VER

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="HMMER"
LABEL software.version="${HMMER_VER}"
LABEL description="Biological sequence analysis using profile hidden Markov models"
LABEL website="http://hmmer.org/"
LABEL license="https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE"
LABEL maintainer="Kutluhan Incekara"
LABEL maintainer.email="kutluhan.incekara@ct.gov"

# get compiled binaries from builder
COPY --from=builder /usr/local/bin/ /usr/local/bin/

ENV LC_ALL=C

# Set the working directory
RUN mkdir -p /opt/work
WORKDIR /opt/work

CMD ["/bin/sh","-c","/bin/bash"]