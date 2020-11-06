# MetaRiboSeq Workflows

This directory contains the workflow for the MetaRiboSeq method.
It uses nextflow as the workflow manager and combination of
Docker and Singularity to build containers.

## Building Images

The `singularity` images for this workflow are built using both `docker`
and `singularity`. Docker is used since it is much more performant and
flexible than `singularity` especially when used with buildkit. Two
docker images are used, one that is a base image intend for broad
use and a second that is specific to metariboseq. The `singularity` image
packages the metariboseq docker image for use on SCG.

See the `./docker` and `./singularity` directories for build scripts etc.

# Nextflow Workflows

`Nextflow` is used to implement the actual workflow; the resulting scripts
can be compared to those developed by `snakemake` to see which is best suited
to our needs going forward. Note that a current version of `nextflow` is used
and that can be found in `/labs/asbhatt/tools/swtools/bin`.

Some common conventions:

1. All nextflow scripts have the suffix `.nf`
2. A shell script, typically `run-nextflow.sh` is used to run `nextflow` with
appropriate parameters. The samples are specified via a parameter file,
as are common workflow options.
3. Nextflow configuration files are specified by these scripts as appropriate
for runing locally or using `slurm`. These are in the `nextflow-configs`
directory.

The parameters file is in `yaml` format and is structured as follows:

```yaml
option1: some options
option2: some options
sampleSpecs:
    - name: <name>
      metagenomic: <metogenomic-data>
      metariboseq: <metariboseq-data>
```

The options are used to configure various command line arguments
passed to the various commands used, e.g. number of threads for `spades`,
or arguments to `trim_galore`. The `sampleSpecs` simple enumerate all
of the samples to be processed by giving each a name and the file name
component (without `.fq.gz`) for each paired set of sequencing files.
Here is a complete example.

```yaml
trimGaloreOptions: "--cores 4 -q 30 --illumina"
spadesOptions: "--threads=4 --memory=96"
alignmentMemory: "96 GB"
assemblyMemory: "96 GB"
bowtieIndexOptions: "--threads 4"
bowtieAlignmentOptions: "--threads 4"
sampleSpecs:
    - name: sampleA
      metagenomic: SampleADNA
      metariboseq: BrayonRibo_1_S7_R1_001
```

All workflows are in the `workflows` directory:

1. `subsample`: generates subsampled input data for testing/development.
2. `metariboseq`: the actual metariboseq workflow.

## Testing and Development

A nextflow script is provided for subsampling a given set of samples to
create a small 'test' sample for development purposes. The full samples
take on the order of 24-48 hours to assemble and hence are cumbersome
to work with. See the `subsample` directory for details. The `run-nextflow.sh`
script runs nextflow with appropriate parameter files and configuration
to generate the subsampled data.

## Metariboseq and Analysis Workflows

The `run-nextflow.sh` scripts accepts one of the following parameters:

- `test`: runs with subsampled test data
- `scg-sampleA`: runs on scg for sampleA
- `scg-all`: runs all samples on scg

Making changes to the parameter files should be fairly self-evident.

## Metariboseq Nextflow Workflow

The `nextflow` execution graph is complicated by the need to create
assemblies for the `metagenomic` files, but not for the `metariboseq`
files and to then build indices and aligments using the `metagenomic`
assembly. The actual `.nf` file is commented, but this is really the
only complication. In terms of understanding the nextflow itself
the documentation and tutorials are (nextflow.io) are reasonably clear.
