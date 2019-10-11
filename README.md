# Bhattlab workflows
Computational workflows for metagenomics tasks, packaged with Snakemake and singularity

### Table of contents

 1. Setup
 2. Running a workflow
 3. Available workflows
	 a. **Preprocessing** metagenomic data
	 b. Metagenomic **Assembly**
	 c. Metagenomic **Binning**
	 f.  Metagenomic classification with **[Kraken2](https://github.com/bhattlab/kraken2_classification)**
	 d. **Sourmash** read comparison
	 e. **Download SRA** data
	 
	 f. Comparative microbial genomics pipelines
	  

### Quickstart


### First time setup


First, install [miniconda3](https://conda.io/miniconda.html) on the cluster, as that's where you should be doing most of your workflow work.

```
# ON SCG
# to download installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# to run installer
sh Miniconda3-latest-Linux-x86_64.sh
# follow instructions, select a place to install miniconda. We recommend installing in the lab directory, /labs/asbhatt/YOURNAME/miniconda3
#install snakemake
conda install -c bioconda -c conda-forge snakemake
```
Then, if you haven't already, set up a Snakemake profile for SCG by following the instructions [here](https://github.com/bhattlab/slurm).
Then clone this github repository to a place on scg. I keep a 'projects' folder in my home directory for cloning repos.
```
# ON SCG
cd ~/projects
git clone https://github.com/bhattlab/bhattlab_workflows.git
cd ~/projects/bhattlab_workflows
```

## Running the workflow
Using the preprocessing workflow here as an example. You'll have to change options in the configuration file to match where your data lives on SCG, etc.
```
snakemake --configfile config_preprocessing.yaml --snakefile path/to/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity --singularity-args '--bind /labs/ --bind /scratch/ '
```

For workflows other than preprocessing, conda environments are used.  Set these up like so:

```
# example for the assembly environment
# where your working directory is the location you cloned the github repository
# ~/projects/bhattlab_workflows  for example
conda env create -f envs/assembly.yaml
```

Then run the workflow as follows:

```
conda activate assembly
snakemake --configfile config_assembly.yaml --snakefile path/to/assembly.snakefile \
--profile scg --jobs 100
```



### Containers
The containers for this repository can be found at https://www.singularity-hub.org/collections/2541. No intervention is required--they are downloaded and used automatically by snakemake when `--use-singularity` is provided as in the example command above.


# Preprocessing

To use this pipeline, edit parameters in the config_preprocessing.yaml, and specify the proper path to config file in the submission script.  The workflow is run with the first snakemake command above.

**The only parameters that need changing in the config file**
1. directory path containing demultiplexed raw fastq files (DATA_DIR)
2. root directory path for output files (PROJECT_DIR)
3. directory path for the host reference genome (BWA index)
4. Directory of the scripts folder from the github, contains sync.py and plot_readcounts.R

*This program runs under the assumption samples are named <sample_id>\_[R]1.fastq[fq].gz and <sample_id>\_[R]2.fastq[fq].gz.* The R1/R2 and suffix must be specified in the config.

**This script will create the following folders:**
- PROJECT_DIR/01_processing/00_qc_reports/pre_fastqc
- PROJECT_DIR/01_processing/00_qc_reports/post_fastqc
- PROJECT_DIR/01_processing/01_trimmed
- PROJECT_DIR/01_processing/02_dereplicate
- PROJECT_DIR/01_processing/03_sync
- PROJECT_DIR/01_processing/04_host_align
- PROJECT_DIR/01_processing/05_sync


The files that can then be used in downstream analyses will be in PROJECT_DIR/01_processing/05_sync/ as {sample}\_1.fq, {sample}\_2.fq, {sample}\_orphans.fq

### Reference genomes for removal of host reads
For Bhatt lab purposes, we only conduct experiments on two hosts, humans and mice. You can specify the host reference genomes in the config using the following directories.
- **Humans:**
``` /labs/asbhatt/data/host_reference_genomes/hg19/hg19.fa ```
- **Mice:**
``` /labs/asbhatt/data/host_reference_genomes/mm10/mm10.fa ```

If working on a different cluster or different model organism, these are the steps necessary to build the host reference index for alignment. I am showing the steps used to build Bhatt lab hosts from above using BWA.

Download reference genomes
```
mkdir -p /labs/asbhatt/data/host_reference_genomes/ # change to preferred directory path
cd /labs/asbhatt/data/host_reference_genomes/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit;
```
Convert to fasta format
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa;
chmod +x twoBitToFa;
./twoBitToFa hg19.2bit hg19.fa;
```
Create bowtie index
```
bwa index hg19.fa
```

# Assembly

Assemble preprocessed read data, using Megahit, Spades or both. Evaluates result with Quast. Accepts a tab-delimited input table, specified in the config.yaml (see example below). Comment lines are specified with the # character. An input table is automatically generated at the end of the preprocessing pipeline and can be found at `01_processing/assembly_input.txt`.

```
#Sample Reads1.fq[.gz][,Reads2.fq[.gz][,orphans.fq[.gz]]]
sample_a    a_1.fq,a_2.fq,a_orphans.fq
sample_b    b_1.fq,b_2.fq,b_orphans.fq
```

# Sourmash
Workflow for kmer-based comparison of sequencing reads. For the details on what sourmash and MinHash are, see the [sourmash docs](https://sourmash.readthedocs.io/en/latest/). Install and activate the `sourmash.yaml` conda environment in the `envs` folder.

There are two inputs required in the config file: the final output directory of the preprocessing pipeline (something like `preprocessing/01_processing/05_sync`) and the output directory. By default, kmer comparisons are calculated at k=21, k=31 and k=51. Output files include:
- Matrices of pairwise jaccard distances between all samples (`04_sourmash_compare/compare_k21.csv`)
- Clustered heatmaps of pairwise distances (`04_sourmash_compare/compare_k21_heatmap_ward.D2.pdf`)

Steps of the workflow include:
- Concatenate all reads into a single file
- Trim lowly abundant kmers
- Build MinHash sketches for each file
- Compare each pair of signatures to get pairwise Jaccard distances
- Downstream processing and plotting

# Metagenomic binning
Binning is essentially clustering for assembled contigs. Create draft metagenome-assembled genomes and evaluate their completeness, contamination and other metrics with these helpful tools!

There are two binning workflows in the `binning` folder. `bin_metabat.snakefile` uses a single binning method ([metabat2](https://peerj.com/articles/1165/)), while `bin_das_tool.snakefile` uses several tools and integrates the result with [DASTool](https://www.nature.com/articles/s41564-018-0171-1). Both use the same downstream evaluation and reporting tools. The DASTool pipeline was made by [Alyssa Benjamin](https://github.com/ambenj).

In contrast to other workflows, you must run binning individually for each sample. Change the following options in the `binning/config_binning.yaml` file to match your project:
- assembly
- sample
- outdir_base
- reads1
- reads2
- read_length

You can launch either workflow using singularity to manage all dependencies with a command like, submitting jobs to the SCG cluster:
```
snakemake --configfile path/to/config_binning.yaml --snakefile path/to/bin_metabat.snakefile \
--profile scg --jobs 100 --use-singularity --singularity-args '--bind /labs/ --bind /scratch/ '
```

## Running binning on many samples
To make this workflow easy to run on many samples, you need to make a configuration file for each one. If all were preprocessed and assembled with our workflows in the same directory, first make one config file to match your samples. Call this `config_example.yaml`. Put your list of samples, one per line, in `sample_list.txt`. Then you can replace the sample name in the configfile in a loop, changing SAMPLE to the name you used with the example configfile:
```
mkdir configfiles
while read line; do
    echo "$line"
    sed "s/SAMPLE/$line/g" config_example.yaml > configfiles/config_"$line".yaml
done < todo_binning.txt 
```

Then, you can run the snakemake workflow for each configfile. Either do this in separate windows in tmux, or run it as a loop. This loop is sequential, but you could even get fancy and run something in parallel with xargs... Change the paths here to correspond to where you have the snakefile and configfiles. 
```
for c in configfiles/*.yaml; do
    echo "starting $c"
    snakemake --snakefile ~/projects/bhattlab_workflows/binning/bin_das_tool.snakefile --configfile "$c" --use-singularity --singularity-args '--bind /labs/ --bind /scratch/' --profile scg --jobs 99 --rerun-incomplete
done
```

# Dumping reads from SRA 
Given a SRR ID number, this will download reads from SRA, using many threads in parallel to speed up the process. Specify the SRR ID number on the command line (replacing SRRXXXXX) . Reads are downloaded into a new folder in your current working directory. 
```
snakemake --snakefile /path/to/sra_download/sra_download.snakefile \
--config srr=SRRXXXXX \
--use-singularity --profile scg --jobs 1
```

# Classification and taxonomic barplots
Deprecated. See our [Kraken2](https://github.com/bhattlab/kraken2_classification) github for the most up to date classification workflow.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE3NTk4MjIwMTddfQ==
-->