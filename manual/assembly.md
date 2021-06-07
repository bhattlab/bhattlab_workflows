## Assembly

Assemble preprocessed shotgun metagenomic data, using Megahit, Spades or both! Results are automatically evaluated with Quast for contiguity. Output of this pipeline easily fits into our [Metagenomic Binning](manual/binning.md) pipeline. Single, paired-end, or paired-end with orphans sequencing data can be used. 

### Setup
To run assembly, first copy the `assembly/config_assembly.yaml` file to a working directory. Edit the config file to chose the assembler you want, the output directory, and the location of a `sample_table` which defines the mapping from sample names to sequence read files. This file must be tab-delimited with two columns-  sample names in the first column, a comma-separated list of paths to sequencing reads in the second column. Comment lines are specified with the # character. An input table is automatically generated at the end of the preprocessing pipeline and can be found at `01_processing/assembly_input.txt`. An example is below.

```
#Sample Reads1.fq[.gz][,Reads2.fq[.gz][,orphans.fq[.gz]]]
sample_a    a_1.fq,a_2.fq,a_orphans.fq
sample_b    b_1.fq,b_2.fq,b_orphans.fq
```

### Run (in the Bhatt lab)
Launch this like any other pipeline, with arguments to submit jobs to the SCG scheduler and use singularity containers.

```
snakemake -s GITHUB/PATH/assembly/assembly.snakefile --configfile config_assembly.yaml --profile scg --jobs 99 --use-singularity --singularity-args '--bind /labs/,/oak/,/home/'
```

If you are satisfied with the assembly pipeline after examining the results, please run the cleanup command to remove unnecessary files and save space on SCG.
```
snakemake cleanup -s GITHUB/PATH/assembly/assembly.snakefile --configfile config_assembly.yaml

```
### Run (not in the Bhatt lab)
You can use this command, but might have to add arguments to submit to your SLURM scheduler, or singularity `bind` arguments as necessary. Configure the number of jobs and cores to fit your available resources. 

```
snakemake -s ~/PATH/TO/assembly.snakefile --configfile config_assembly.yaml --use-singularity --jobs 8 --cores 8
```
