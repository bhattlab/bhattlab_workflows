
## Metagenomic binning
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

### Running binning on many samples
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

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTExMjE5NjEyN119
-->