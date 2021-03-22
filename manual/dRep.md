# dRep: deReplication of MAG sets

dRep is a tool developed by Matt Olm. See the [Documentation](https://drep.readthedocs.io/en/latest/) or [paper](https://www.nature.com/articles/ismej2017126) for more info. 

dRep runs pairwise comparisons between all genomes in your input, creates trees, and produces a set of unique (de-replicated) genomes in the final output. This pipeline is designed to work directly from the output of our DAS_Tool binning pipeline. Only a few variables need to specified in the configfile, including the desired output directory and the directory with results from DAS_Tool. 

Run this pipeline the same as any other, specifying the snakefile and configfile on the command line. Add `--use-singularity` to use the singularity image for dRep if you don't have it installed. 

If you want to run dRep on a set of genomes that didn't come from DAS_Tool, it's probably easier to just run the program by itself rather than make your data fit this snakefile. You can use the commands from the dRep rule in the snakefile as an example. 
