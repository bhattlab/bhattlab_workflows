#!/bin/bash

root="$(git rev-parse --show-toplevel)/mge-hotspots"
params="${root}/workflows/params"

#container="${root}/singularity/$1"
#configs="-c ${root}/nextflow-configs/singularity.config"
#containerFlags=-with-singularity "${container}"

workdir=/labs/asbhatt/cosn/nextflow-workdir
paramsfile=""
case "$1" in
"mdurrant")
    paramsfile="${params}/params-mdurrant.yml"
    ;;
"mdurrant-cos-local")
    paramsfile="${params}/params-mdurrant-cos-local.yml"
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
    echo "unrecognised mode: $1, use one of mdurrant, mdurrant-cos-local, stm20x, stm20x-cos-local"
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
