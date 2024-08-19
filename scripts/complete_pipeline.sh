#!/usr/bin/env bash

# Defining command line arguments
while getopts i:o:s: flag
do
    case "${flag}" in
        i) input_path=${OPTARG};;
        o) output_path=${OPTARG};;
        s) path_to_scripts=${OPTARG};;
    esac
done

# /Volumes/Data/PhenoImager_cores/SCAN-B_TNBC_TMA_1A

# Getting input files from input path
input_files=$(ls ${input_path}/*component_data.tif | head -1 | tr "\n" ":")

# Generating input directory and subdirectories
mkdir ${output_path}

mkdir ${output_path}/segmentation_masks
mkdir ${output_path}/original_phenotypes
mkdir ${output_path}/refined_phenotypes

# Running cell segmentation and preliminary cell type assignment
python ${path_to_scripts}/segmentation_and_phenotyping.py -p ${input_files} -o ${output_path}/original_phenotypes -s ${output_path}/segmentation_masks

# Running segmentation refinement

# Getting paths to original phenotypes
original_phenotypes=$(ls ${output_path}/original_phenotypes/*intensity.csv)

# Iterate through files and perform refinement
for element in ${original_phenotypes[@]};
    do Rscript ${path_to_scripts}/refining_phenotypes_GA.R -i ${element} -p ${output_path}/refined_phenotypes