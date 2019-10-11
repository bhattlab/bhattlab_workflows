
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

 - `--profile scg` tells snakemake to submit batch jobs for compute-heavy tasks ([install the SCG profile](https://github.com/bhattlab/slurm)).
 - `--jobs 100` will submit up to 100 jobs at the same time
 - `--use-singularity` tells snakemake to use singularity images for software, if they're specified in the workflow
 - `--singularity-args '--bind /labs/  --bind /home/'` lets singularity access other directories on the cluster

```
snakemake --configfile config_preprocessing.yaml \
--snakefile ~/projects/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--profile scg --jobs 100 --use-singularity \
--singularity-args '--bind /labs/  --bind /home/'
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIyMDk0MjkwNF19
-->