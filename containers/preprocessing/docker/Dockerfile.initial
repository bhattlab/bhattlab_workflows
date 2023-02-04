# syntax=docker/dockerfile:1
FROM ubuntu:20.04

ENV LC_ALL=C
RUN apt-get update -y
RUN apt-get install -y curl unzip
RUN apt-get clean
RUN apt-get autoclean
RUN apt-get install -y ca-certificates

# fix this time zone stuff to install only the tzdata package for
# PST.
ENV TZ=US/Pacific
RUN apt-get install -y add-apt-key
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN useradd -m -d /opt/labmember -s /bin/bash labmember
USER labmember
WORKDIR /opt/labmember
RUN mkdir bin
ENV HOME=/opt/labmember

# install miniconda and then bioconda.
RUN curl -so Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/labmember/miniconda3 
RUN rm -f ./Miniconda3-latest-Linux-x86_64.sh

ENV PATH=$PATH:$HOME/bin:$HOME/miniconda3/bin
RUN <<EOT
conda update -y -n base -c defaults conda
conda init bash
EOT

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

RUN conda env create -f preprocessing.yml

RUN conda clean -afy
