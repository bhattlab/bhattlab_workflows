#!/bin/bash

root="/labs/asbhatt/cosn/bhattlab_workflows/metariboseq"

container="${root}/singularity/$1"
nextflow=$2
shift
shift
params="${root}/workflows/params"

configs="-c ${root}/nextflow-configs/singularity.config"
paramsfile=""
case "$1" in
"test")
    paramsfile="${params}/params-test-subsampled.yml"
    ;;
"scg-sampleA")
    configs="$configs,${root}/nextflow-configs/slurm.config"
    paramsfile="${params}/params-sampleA.yml"
    ;;
"scg-all")
    configs="$configs,${root}/nextflow-configs/slurm.config"
    paramsfile="${params}/params-all.yml"
    ;;
*)
    echo "unrecognised mode: $1"
    exit 1
    ;;
esac

sd=./reports
mkdir -p $sd
td="$(date +"%Y_%m_%d:%H:%M")"
debugOptions="-dump-hashes -with-trace $sd/trace-${td}.txt -with-report $sd/report-${td}.html -with-timeline $sd/timeline-${td}.html"
debugOptions="-dump-hashes"
debugOptions=""

set -x
nextflow run -work-dir /labs/asbhatt/cosn/nextflow-workdir \
    $debugOptions \
    $configs \
    -resume \
    -params-file ${paramsfile} \
    -with-singularity ${container} \
    ${nextflow}
