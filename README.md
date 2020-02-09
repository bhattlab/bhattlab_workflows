# Bhattlab workflows
Computational workflows for metagenomics tasks, packaged with Snakemake and singularity

### Table of contents

 1. [Setup](manual/setup.md)
 2. [Running a workflow](manual/running.md)
 3. Available workflows
    - [**Preprocessing** metagenomic data](manual/preprocessing.md)
    - [Metagenomic **Assembly**](manual/assembly.md)
    - [Metagenomic **Binning**](manual/binning.md)
    - [Metagenomic classification with **Kraken2**](https://github.com/bhattlab/kraken2_classification)
    - [**Sourmash** read comparison](manual/sourmash.md)
    - [**Download SRA** data](manual/download_sra.md)
    - [Comparative microbial genomics pipelines](manual/comparative_genomics.md)
	  

### Quickstart
```
snakemake --configfile config_preprocessing.yaml \
--snakefile ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity \
--singularity-args '--bind /labs/,/oak/,/home/'
```
