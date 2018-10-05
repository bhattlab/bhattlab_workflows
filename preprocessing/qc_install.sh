#!/usr/bin/env bash

IS_CONDA="$(which conda)"

if [ ${#IS_CONDA} -eq "0" ]
    then
        echo "It seems that you do not have conda installed. Download it from https://www.continuum.io/downloads"
        exit 1
fi

echo "Removing preexisting qc environment if exists"
source deactivate
conda remove --name qc_sm --all --yes;

echo "Creating the new qc environment from qc_env.yaml"
conda env create -f qc_env.yaml;

if [ $? -ne 0 ]
    then
        echo "An issue occured creating the conda environment. Please consider manually creating a python virtualenv and updating proper packages with the requirements.txt file."
        exit 1
fi

################################################################################
echo "------------------------------------------------------------------------------------------------------"
echo "INSTALLATION COMPLETE"
echo "------------------------------------------------------------------------------------------------------"
echo "NOTE: You must specify the required information in the config.yaml file to run the workflow properly."
echo "NOTE: Once configured, enter:"
echo ""
echo "> source activate qc_sm"
echo "> snakemake"
echo ""
echo "In this directory to run the workflow."
echo "------------------------------------------------------------------------------------------------------"
