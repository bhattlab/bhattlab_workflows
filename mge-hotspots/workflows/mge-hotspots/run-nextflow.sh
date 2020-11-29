#!/bin/bash

root="$(git rev-parse --show-toplevel)/mge-hotspots"
params="${root}/workflows/params"

#container="${root}/singularity/$1"
#configs="-c ${root}/nextflow-configs/singularity.config"
#containerFlags=-with-singularity "${container}"

workdir=/labs/asbhatt/cosn/nextflow-workdir
paramsfile=""
case "$1" in
"rare-insertions")
    paramsfile="${params}/params-rare-insertions.yml"
    ;;
"rare-insertions-cos-local")
    paramsfile="${params}/params-rare-insertions-cos-local.yml"
    workdir=./workdir
    ;;
"stm20x")
    paramsfile="${params}/params-stm_20x.yml"
    ;;
"stm20x-cos-local")
    paramsfile="${params}/params-stm_20x-cos-local.yml"
    workdir=./workdir
    ;;
*)
    echo "unrecognised mode: $1, use one of 'rare-insertions' or 'salmonella'"
    exit 1
    ;;
esac

set -x
nextflow run -work-dir "${workdir}" \
    $configs \
    -resume \
    -params-file "${paramsfile}" \
    ${containerFlags} \
   mge-hotspots.nf
