## Metagenomic binning
Binning is essentially clustering for assembled contigs. Create draft metagenome-assembled genomes and evaluate their completeness, contamination and other metrics with these helpful tools!

This pipeline uses [DASTool](https://www.nature.com/articles/s41564-018-0171-1), which aggregates results from several binning methods to produce a higher-quality bin set. The DASTool pipeline here uses three binners under the hood: Metabat2, Maxbin and CONCOCT. While we chose these binners and DASTool, there are several other binners and aggregation tools out there in the literature.

Kraken2 is currently used to make a taxonomic assignment of a bin. You need a functional Kraken2 database for this pipeline to work. 

### Running DASTool
This pipeline will work for a set of many samples at once. Copy the `binning/config_binning_manysamp.yaml` to your working directory and change the parameters. 

**Specify the following fields in the configfile (defaults are provided)**
    - outdir_base (a folder for each sample will be created within here)
    - sample_file (see below)
    - read_length 
    - kraken2db
    - custom_taxonomy
    - long_read

Your sample file must be a tab delimited file with three columns and one row for each sample. Forward and reverse reads should be in the same column separated by commas.
```
SAMPLE_NAME    ASSEMBLY   READS1,READS2
```

Then, you can launch the workflow to submit jobs to the SCG SLURM cluster. Note the use of `--use-conda` on the end to get certain rules to work correctly. 
```
snakemake --configfile config_binning_manysamp.yaml --snakefile path/to/bin_das_tool_manysamp.snakefile \
--profile scg --jobs 999 --use-singularity --singularity-args '--bind /oak/,/labs/,/home' --use-conda
```

### Comparing bins against references
**THIS PIPELINE HAS BEEN DEPRECATED WITHIN THE BHATT LAB, I DELETED THE DATABASES**
To compare your freshly minted bins against reference genomes from Genbank and MAG databases, run the `compare_bins_references.snakefile` pipeline. For each sample, this will compare the bin against the genbank and MAG databases with MASH. The top 50 matches will be passed into fastani, and the top 10 from there will get more careful nucmer alignments and dotplots. This pipeline relies on the directory structure output from the binning pipeline. The file `config_compare_bins_references.yaml` takes two inputs: 
 - base output directory from binnnig. expects one folder for each sample within this directory. This should be where you ran DASTool or Metabat
 - sample file: one column, newline delimited for each sample this is just the sample name. There should be a folder in outdir_base for each sample name specified here. 

### Running on Nanopore assemblies
Change the `long_read` specification in the configfile and put your nanopore reads in the third column (just once). If you also have high-quality short read data, you could run this pipeline in the standard mode by providing your nanopore assembly and the short read pairs. This might be more accurate for depth/coverage calculations, but I'm not sure. 

*Still in development, please file issues on GitHub if you encounter any.*
