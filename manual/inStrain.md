# inStrain 

inStrain is a tool developed by Matt Olm. See the [Documentation](https://instrain.readthedocs.io/en/latest/) or [paper](http://www.nature.com/articles/s41587-020-00797-0) for more info. inStrain has two main uses, which correspond to the two steps in the program. 

## inStrain profile
This step _profiles_ sequencing reads from one sample against a reference genome. It calls SNPs and reports metrics on coverage depth and breadth. If you're just looking for stats on how well a new sample matches a reference (or assembly/MAG), you can just run this step. 

## instrain compare
This step _compares_ multiple runs of profile from different samples against the same reference genome in a pairwise manner. It compares SNPs called in different samples to generate the consensus and population ANI metrics. If you're interested in comparing strains in multiple samples, run this step as well. 

## Running the pipeline
There are two ways to run this pipeline, depending on your data format and research question. 

 1. You have the output of dRep, containing a set of unique MAGs from your sample collection, and you want to compare multiple short-read samples using all of your MAGs as reference. 
 2. You have a SINGLE reference genome (perhaps something from RefSeq, or a MAG) and you want to compare multiple short-read samples using this reference 

### Version 1: output of dRep into inStrain
In the configfile, provide a file for _sample_reads_, the desired output directory, and the directory with the results from dRep. You can provide a _sample_groups_ file if you want possible transmission events to be annotated in the results (groups being different patients, for example). 

Run the `instrain_profile.snakefile` pipeline followed by the `instrain_compare.snakefile` pipeline, using the same configfile. Note the use of the `--use-conda` argument here.

```
# add other arguments as needed, like: --profile scg --jobs 999 -k --singularity-args '--bind /labs/,/oak/,/home/' 
snakemake -s ~/projects/bhattlab_workflows/instrain/instrain_profile.snakefile --configfile config_instrain.yaml --use-singularity --use-conda 

# wait a long time for profile to complete
snakemake -s ~/projects/bhattlab_workflows/instrain/instrain_compare.snakefile --configfile config_instrain.yaml --use-singularity --use-conda 
```

User-defined settings can include: 
 1. instrain_min_reads: Minimum number of reads to consider a sample in profile and compare steps. The default is set at 20,000 reads, which eliminates low-coverage comparisons and saves a lot of time. This is appropriate for for bacterial genomes, but should be changed if you have smaller genomes (phage for example). 
 2. instrain_max_reads: Above this number of reads, bam files will be randomly subsampled to save time when coverage is extremely high. 
 3. cluster_file: By default, genomes from dRep that have >1x coverage in >=2 samples are considered for inStrain profiling. You can override this by providing a file here with names and paths to the genomes. 

Running inStrain can be very time intensive, especially when the number of clusters is high (>1000). This is because each sample has to be profiled against each cluster (n\*m time), and then pairwise comparisons between all samples have to be done for each cluster (n\*n\*m time!) Try to compare against a smaller number of clusters, if possible. You can also provide a set of integers to the _limit_clusters_ argument in the configfile like "0,1,2,3" which will limit the comparison to the specified clusters. You could make a configfile for the first 100 clusters, and so on, in this fashion. 


### Version 2: profile and compare samples against a single MAG 
Run this the same way as above, but provide a file to the config option "cluster_file". This file should contain the name and path of the genome you want to map reads against. You can ignore the config arguments about dRep as they won't be used. If you want to compare _multiple_ reference genomes that didn't come from dRep, I would recommend setting this up several times, one for each reference. 

Note the use of the `--use-conda` argument here.
```
# add other arguments as needed, like: --profile scg --jobs 999 -k --singularity-args '--bind /labs/,/oak/,/home/' 
snakemake -s ~/projects/bhattlab_workflows/instrain/instrain_profile_single.snakefile --configfile config_instrain.yaml --use-singularity --use-conda 

# wait a long time for profile to complete
snakemake -s ~/projects/bhattlab_workflows/instrain/instrain_compare_single.snakefile --configfile config_instrain.yaml --use-singularity --use-conda 
```

For example, If you generate an assembly of a bloodstream isolate from sequencing reads, and you want to compare reads from many stool/metagenomic samples to the isolate, you should ALSO include the isolate sequencing reads in your comparison. This will allow you to get metrics on strain diversity in the isolate, as well as make the pairwise comparisons between the isolate and the stool sample very easy. 