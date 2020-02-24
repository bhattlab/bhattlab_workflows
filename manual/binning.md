## Metagenomic binning
Binning is essentially clustering for assembled contigs. Create draft metagenome-assembled genomes and evaluate their completeness, contamination and other metrics with these helpful tools!

The current best practice for binning is to use [DASTool](https://www.nature.com/articles/s41564-018-0171-1), which aggregates results from several binning methods to produce a higher-quality bin set. The DASTool pipeline here uses three binners at its heart: Metabat2, Maxbin and CONCOCT. 

### Running DASTool
DASTool can now be run across many samples at once! Use `binning/bin_das_tool_manysamp.snakefile` and specify parameters in `binning/config_binning_manysamp.yaml`

**Specify the following fields in the configfile (defaults are provided)**
    - outdir_base (a folder for each sample will be created within here)
    - sample_file (see below)
    - read_length 
    - kraken2db
    - custom_taxonomy
    - long_read

Your sample file must be a tab delimited file with three columns and one row for each sample. Each row should look like this, with tab as the delimiter. Assembly and reads are file paths, and forward/reverse reads are specified with commas.
```
SAMPLE_NAME    ASSEMBLY   READS1,READS2
```

Then, you can launch the workflow to submit jobs to the SCG SLURM cluster.
```
snakemake --configfile config_binning_manysamp.yaml --snakefile path/to/bin_das_tool_manysamp.snakefile \
--profile scg --jobs 999 --use-singularity --singularity-args '--bind /oak/,/labs/,/home'
```

### Running on Nanopore assemblies
Change the `long_read` specification in the configfile and put your nanopore reads in the third column (just once). If you also have high-quality short read data, you could run this pipeline in the standard mode by providing your nanopore assembly and the short read pairs. This might be more accurate for depth/coverage calculations, but I'm not sure. 

*Still in development, please file issues on GitHub if you encounter any.*