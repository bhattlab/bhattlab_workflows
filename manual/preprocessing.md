
## Preprocessing

Preprocessing raw metagenomic data is necessary before any other application. Our preprocessing pipeline does the following tasks:

 0. Initial quality control check 
 1. Trimming of low quality bases and discarding of reads that fall below a length limit
 2. Deduplicaiton of exactly matching sequencing reads
 3. Removal of reads that align against the host genome
 4. Final quality control check

To use this pipeline, copy `preprocessing/config_preprocessing.yaml` to your working directory and edit the options to match your project. The workflow is run with the first snakemake command below.
```
snakemake --configfile config_preprocessing.yaml \
--snakefile ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity \
--singularity-args '--bind /labs/  --bind /home/'
```

### Change configfile options
1. directory path containing demultiplexed raw fastq files (raw_reads_directory)
2. root directory path for output files (output_directory)
3. How paired end reads are specified (read_specification)
4. Read file extension (extension)
5. Path for the host reference genome (host_genome)


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
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEyODI3MTExNzYsNjQ4MTcwODA2XX0=
-->