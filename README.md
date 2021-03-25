# Bhattlab workflows
Computational workflows for metagenomics tasks, packaged with Snakemake, Singularity and Conda.

### Table of contents

 1. [Setup](manual/setup.md)
 2. [Running a workflow](manual/running.md)
 3. Available workflows
    - [**Preprocessing** metagenomic data](manual/preprocessing.md)
    - [Metagenomic **Assembly**](manual/assembly.md)
    - [Metagenomic **Binning**](manual/binning.md)
    - [**DeReplication** of binned genomes](manual/dRep.md)
    - [**inStrain** strain-diversity aware comparison of samples](manual/inStrain.md)
    - [Metagenomic classification with **Kraken2**](https://github.com/bhattlab/kraken2_classification)
    - [**Sourmash** read comparison](manual/sourmash.md)
    - [**Download SRA** data](manual/download_sra.md)
	- [**ARG detection** with RGI](manual/arg.md) 
    - [**Viral** contig prediction](manual/viral.md) 
    - [Comparative microbial genomics pipelines](manual/comparative_genomics.md)
    - [MetaRibo-Seq](metariboseq/README.md)

### Quickstart
If you're in the Bhatt lab and working on SCG, this command is an example of how to run the workflows. Other users will need to change these options (see [Running a workflow](manual/running.md))
```
snakemake --configfile config_preprocessing.yaml \
--snakefile ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity \
--singularity-args '--bind /labs/,/oak/,/home/'
```

**Important** after running the preprocessing or assembly pipelines, run the cleanup rule in the snakefile. This will delete all the unnecessary files and save us space on SCG! 
```
snakemake cleanup -s PATH/TO/SNAKEFILE --configfile YOUR_CONFIGFILE.yaml

```

You can also remove the `.snakemake` folder, which only contains metadata and cache from the run. 
```
rm -rf .snakemake 
```
