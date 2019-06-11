## First time setup

First, install [miniconda3](https://conda.io/miniconda.html) on the cluster, as that's where you should be doing most of your workflow work.

```
# ON SCG
# to download installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# to run installer
bash https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# follow instructions, select a place to install miniconda. We recommend installing in the lab directory, /labs/asbhatt/YOURNAME/miniconda3
```
Then, if you haven't already, set up a Snakemake profile for SCG by following the instructions [here](https://github.com/bhattlab/slurm). 
Then clone this github repository to a place on scg. I keep a 'projects' folder in my home directory for cloning repos.
```
# ON SCG
cd ~/projects
git clone git@github.com:bhattlab/bhattlab_workflows.git
cd ~/projects/bhattlab_workflows
```
Then create an environment for your workflow of choice. For example, to start with the preprocessing workflow:

```
conda env create -f envs/preprocessing.yaml
```

## Running the workflow
Using the preprocessing workflow here as an example. You'll have to change options in the configuration file to match where your data lives on SCG, etc.
```
source activate preprocessing
snakemake --configfile config_preprocessing.yaml --snakefile preprocessing.snakefile --profile scg --jobs 10
```


# Preprocessing 

To use this pipeline, edit parameters in the config_preprocessing.yaml, and specify the proper path to config file in the submission script.

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

You can launch either workflow using singularity images with a command like 

# Classification and taxonomic barplots
Deprecated. See our [Kraken2](https://github.com/bhattlab/kraken2_classification) github for the most up to date classification workflow.
