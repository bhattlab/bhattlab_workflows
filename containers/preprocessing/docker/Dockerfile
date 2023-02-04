# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:focal

USER root
RUN apt-get update && \
  apt-get install -yq tzdata && \
  ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime && \
  dpkg-reconfigure -f noninteractive tzdata
USER $MAMBA_USER

RUN <<EOF cat > preprocessing.yml
name: preprocessing
channels:
  - bioconda
  - conda-forge
dependencies:
  - ncurses
  - python=3.10
  - fastqc=0.11.9
  - multiqc=1.13
  - htstream=1.3.3
  - bwa=0.7.17
  - samtools=1.16.1
  - trim-galore=0.6.7
EOF

RUN micromamba install -y -n base -f preprocessing.yml && \
  micromamba clean --all --yes