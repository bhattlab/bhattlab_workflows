FROM ghcr.io/cosnicolaou/bhattlab-ubuntu-focal-base

# install linux packages
USER root
RUN apt-get update
RUN apt-get install -y fastqc=0.11.9+dfsg-2
RUN apt-get install -y samtools=1.10-3
RUN apt-get install -y bamtools=2.5.1+dfsg-5build1
RUN apt-get clean
RUN apt-get autoclean

# install bioconda packages
USER labmember
WORKDIR /opt/labmember
ENV HOME=/opt/labmember
RUN PATH=$HOME/miniconda3/bin:$PATH; \
	conda create -n metariboseq -y cutadapt && \
	conda clean -a -y

# install packages with custom installers and configuration

# trim_galore
RUN curl -fsSL -o trim_galore.tar.gz \
	https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz && \
	tar xvzf trim_galore.tar.gz && \
	mv TrimGalore-0.6.6/trim_galore bin && \
	rm -rf ./TrimGalore-0.6.6 trim_galore.tar.gz

# bowtie 1 and 2.
RUN curl -fsSL -o bowtie.zip \
	https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip/download && \
	unzip bowtie.zip && \
	rm bowtie.zip
RUN curl -fsSL -o bowtie2.zip \
	https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip/download && \
	unzip bowtie2.zip && \
	rm  bowtie2.zip

# spades assembler
RUN curl -fsSL -o spades.tgz \
	http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz && \
	tar zxf spades.tgz && \
	rm spades.tgz

ENV PATH=$HOME/SPAdes-3.14.1-Linux/bin:$HOME/bowtie-1.3.0-linux-x86_64:$HOME/bowtie2-2.4.2-linux-x86_64:$PATH
