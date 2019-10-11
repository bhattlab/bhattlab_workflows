
## First time setup
If you're in the Bhatt lab, most of your work will take place on the SCG cluster. Our [scg_tools](https://github.com/bhattlab/scg_tools) repository has all the information you need to get set up with cluster computing. 

If you're not in the Bhatt lab, you can still use all the workflows here. The file paths to databases and such will need to be changed, but most things are packaged with Singularity which means you won't need to install any custom software. 

### Installing miniconda3 and Snakemake
First, install [miniconda3](https://conda.io/miniconda.html) on the cluster, as that's where you should be doing most of your workflow work.
```
# ON SCG
# to download installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# to run installer
sh Miniconda3-latest-Linux-x86_64.sh
# follow instructions, select a place to install miniconda. We recommend installing in the lab directory, /labs/asbhatt/YOURNAME/miniconda3
#install snakemake
conda install -c bioconda -c conda-forge snakemake
```
Then, if you haven't already, set up a Snakemake profile for SCG by following the instructions [here](https://github.com/bhattlab/slurm).
Then clone this github repository to a place on scg. I keep a 'projects' folder in my home directory for cloning repos.
```
# ON SCG
mkdir -p ~/projects && cd ~/projects
git clone https://github.com/bhattlab/bhattlab_workflows.git
cd bhattlab_workflows
```
That's it! Now you're ready to start [running the workflows.](manual/running.md)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNzEzMjEzOTM0XX0=
-->