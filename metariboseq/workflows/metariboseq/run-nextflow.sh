#!/bin/bash

configs="-c ../nextflow-configs/singularity.config"
paramsfile=""
case "$1" in
"test")
    paramsfile="params-test-subsampled.yml"
    ;;
"scg-sampleA")
    configs="$configs,../nextflow-configs/slurm.config"
    paramsfile="params-sampleA.yml"
    ;;
"scg-all")
    configs="$configs,../nextflow-configs/slurm.config"
    paramsfile="params-all.yml"
    ;;
*)
    echo "unrecognised mode"
    exit 1
    ;;
esac

sd=./reports
mkdir -p $sd
td="$(date +"%Y_%m_%d:%H:%M")"
debugOptions="-dump-hashes -with-trace $sd/trace-${td}.txt -with-report $sd/report-${td}.html -with-timeline $sd/timeline-${td}.html"
debugOptions="-dump-hashes"
debugOptions=""
nextflow run -work-dir /labs/asbhatt/cosn/nextflow-workdir $debugOptions $configs -resume -params-file ${paramsfile} metariboseq.nf
