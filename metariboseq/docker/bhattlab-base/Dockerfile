FROM ubuntu:20.04

ENV TZ=US/Pacific
RUN apt-get update -y
RUN apt-get install -y curl tzdata unzip pigz
RUN apt-get clean
RUN apt-get autoclean
RUN useradd -m -d /opt/labmember -s /bin/bash labmember
USER labmember
WORKDIR /opt/labmember
RUN mkdir bin

# install miniconda and then bioconda.
RUN curl -so Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	sh ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/labmember/miniconda3 && \
	rm -f ./Miniconda3-latest-Linux-x86_64.sh

ENV HOME=/opt/labmember
RUN PATH=$HOME/miniconda3/bin:$PATH; \
	conda config --add channels defaults && \
	conda config --add channels bioconda && \
	conda update -y -n base -c defaults conda && \
	conda init bash
