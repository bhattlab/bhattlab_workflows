
## Dumping reads from SRA 
Given a list off SRR IDs, this uses fasters-dump to download reads from SRA, using many threads in parallel to speed up the process. Specify the list in the config file. Reads are downloaded into the current working directory.
```
snakemake -s /path/to/sra_download/sra_download.snakefile
 --configfile config_sra.yaml --cores 1 --profile scg —jobs 1
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTkzMTI1ODQ1Ml19
-->
