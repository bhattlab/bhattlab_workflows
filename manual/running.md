## Running a workflow
Most workflows have two components: 

 - A *snakefile* which specifies the commands to run in the workflow
 - A *configfile* which controls options, input/output directories, databases, sample names, etc.

Running a workflow can be as easy as modifying the *configfile* and typing a command!

### Example
We'll use the [preprocessing](preprocessing.md) workflow here as an example. Let's assume I have some new metagenomic sequencing data stored in `/labs/ashbatt/bsiranos/exciting_new_project/raw_data`.

I'll work within the `exciting_new_project` directory for the preprocessing pipeline. First, go to that directory and copy over the example configfile from where you cloned the github repository (in [setup](setup.md)). Then edit the configfile to specify the right paths and options.
```
cd  /labs/ashbatt/bsiranos/exciting_new_project/
cp ~/projects/bhattlab_workflows/preprocessing/config_preprocessing.yaml .
vim config_preprocessing.yaml
```
Run snakemake, giving it the path to the workflow snakefile, the configfile, and some other options:

 - `--profile scg` tells snakemake to submit batch jobs for compute-heavy tasks ([install the SCG profile](https://github.com/bhattlab/slurm)). Only use this if you are in the Bhatt lab and working on the SCG cluster. 
 - `--jobs 100` will submit up to 100 jobs at the same time. If you're not using the profile option above, change this to the number of cores you want to use on the machine, and add a `--cores X` option as well. So `--jobs 16 --cores 16` will run up to 16 jobs simultaneously using up to 16 cores.
 - `--use-singularity` tells snakemake to use singularity images for software, if they're specified in the workflow
 - `--use-conda` may be necessary, check the workflow specific page. 
 - `--singularity-args '--bind /labs/,/oak/,/home/'` lets singularity access other directories on the cluster. This is only necessary for running on SCG in the Bhatt lab. Other users might need to change these directories 

```
snakemake --configfile config_preprocessing.yaml \
--snakefile ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity \
--singularity-args '--bind /labs/,/oak/,/home/'
```
