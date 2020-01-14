## Preprocessing
**NEW 2020-01-14**
The preprocessing pipeline has been reworked significantly! _hts\_SuperDeduper_ is used to do deduplication in one pass, which happens before trimming. And now that I know about the `-s` argument to `samtools fastq`, there's no longer _any_ syncing needed! Benchmarking has also been added to each step.

Preprocessing raw metagenomic data is necessary before any other application. Our preprocessing pipeline does the following tasks:

 0. Initial quality control check 
 1. Deduplicaiton of exactly matching sequencing reads - PCR artifacts
 2. Trimming of low quality bases and discarding of reads that fall below a length limit
 3. Removal of reads that align against the host genome
 4. Final quality control check

To use this pipeline, copy `preprocessing/config_preprocessing.yaml` to your working directory and edit the options to match your project. Assuming you have this github cloned into a directory called `~/projects`, The workflow can be run with a snakemake like this. You must set up the (scg cookiecutter)[https://github.com/bhattlab/slurm] profile to use the `--profile scg` option to submit batch jobs to the SLURM cluster. 
```
snakemake -s ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--configfile config_preprocessing.yaml \
--profile scg \
--jobs 100 \
--use-singularity \
--singularity-args '--bind /labs,/oak,/home'
```
**Important** after running the pipeline, run the cleanup rule in the snakefile. This will delete all the unnecessary fastq files and save us space on SCG! 
```
snakemake cleanup \
-s ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--configfile config_preprocessing.yaml

```

### Change configfile options
1. directory path containing demultiplexed raw fastq files (raw_reads_directory)
2. root directory path for output files (output_directory)
3. How paired end reads are specified (read_specification)
4. Read file extension (extension)
5. Path for the host reference genome (host_genome)
6. Quality trimming settings (should be fine as is)


*This program runs under the assumption samples are named <sample_id>\_[R]1.fastq[fq].gz and <sample_id>\_[R]2.fastq[fq].gz.* The R1/R2 and suffix must be specified in the config.

**This script will create the following folders:**
- PROJECT_DIR/01_processing/00_qc_reports/pre_fastqc
- PROJECT_DIR/01_processing/00_qc_reports/post_fastqc
- PROJECT_DIR/01_processing/01_dedup
- PROJECT_DIR/01_processing/02_trimmed
- PROJECT_DIR/01_processing/05_sync
_The 03_sync and 04_host_align folder have been removed because of improvements in the pipeline, but the numbering convention has not changed to keep consistency with older versions_ 

The files that can then be used in downstream analyses will be in PROJECT_DIR/01_processing/05_sync/ as {sample}\_1.fq.gz, {sample}\_2.fq.gz, {sample}\_orphans.fq.gz

### Examine your output!
Check to make sure your data is high quality after running this pipeline! View the multiqc reports for a good overview of quality. Are there any problems indicated? Did adapter sequence contamination get removed from pre to post? Also check the `01_processing/readcounts.pdf` and `01_processing/readcounts.tsv` file for the number of sequencing reads that were kept/discarded at each step. This is helpful for diagnosing problems with data.

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
