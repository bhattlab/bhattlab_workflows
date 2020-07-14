
## Dumping reads from SRA 
Given a list off SRR IDs, this uses fasterq-dump to download reads from SRA, using many threads in parallel to speed up the process. Specify a list of SRA IDs in a newline delimited file, and specify this file as a config argument like in the example below. By default 4 threads are used per download, so change the `cores` and `jobs` parameters to the number of concurrent downloads you want to run at the same time. You can also submit these as batch jobs (add `--profile scg`) if you want to run many in parallel. 

```
# this will run one job with 4 threads at a time
snakemake -s /path/to/sra_download/sra_download.snakefile \
--config srr_list=srr_ids.txt --jobs 4 --cores 4

# srr_ids.txt looks like this: 
SRR12185738
SRR12185739
SRR12185740
```

Note that there is an error that crops up occasionally where a valid SRA ID is supposedly invalid. The only fix right now is to try again, see https://github.com/ncbi/sra-tools/issues/215