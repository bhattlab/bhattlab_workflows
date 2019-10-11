
## Dumping reads from SRA 
Given a SRR ID number, this will download reads from SRA, using many threads in parallel to speed up the process. Specify the SRR ID number on the command line (replacing SRRXXXXX) . Reads are downloaded into a new folder in your current working directory. 
```
snakemake --snakefile /path/to/sra_download/sra_download.snakefile \
--config srr=SRRXXXXX \
--use-singularity --profile scg --jobs 1
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTkzMTI1ODQ1Ml19
-->