# Classification and barplots

Given one or more input datasets, this workflow performs taxonomic classification with kraken, then visualizes 
the resulting compositions in a barplot.

## First time setup

First, install [miniconda3](https://conda.io/miniconda.html).  Then, if you haven't already, set up a Snakemake 
profile for SCG by following the instructions [here](https://github.com/bhattlab/slurm).  Then:

```
conda env create -f /path/to/classification.yaml
```

## Running the workflow

```
source activate classification
snakemake --configfile /path/to/config.yaml -s /path/to/classification/Snakefile --profile scg
```

