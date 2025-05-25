#!/bin/bash

source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

mamba env create -f environment.yml || true

mamba activate alzheimer-analysis

mkdir -p ./data
mkdir -p ./outputs
mkdir -p ./results

wget -O ./data/alzheimer_mouse_data.tar.gz "https://osf.io/download/6833661cfe08c0f335d6a716/"
tar -xzf ./data/alzheimer_mouse_data.tar.gz -C ./data --strip-components=1

# Remove the tar.gz file to save space
rm ./data/alzheimer_mouse_data.tar.gz

python run_analysis.py