#!/bin/bash
set -e

# Get the absolute path of the script directory
SCRIPT_DIR="/app"

# Create output directories
mkdir -p "${SCRIPT_DIR}/data"
mkdir -p "${SCRIPT_DIR}/outputs"
mkdir -p "${SCRIPT_DIR}/results"

# Download data if it doesn't exist
echo "Downloading comparative genomics data..."
if [ ! -f "${SCRIPT_DIR}/data/micrococcus_data.tar.gz" ]; then
    wget -O "${SCRIPT_DIR}/data/micrococcus_data.tar.gz" "https://osf.io/download/nz3sq/"
fi

# Extract data if not already extracted
if [ ! -d "${SCRIPT_DIR}/data/genome_1" ]; then
    echo "Extracting data..."
    tar -xzf "${SCRIPT_DIR}/data/micrococcus_data.tar.gz" -C "${SCRIPT_DIR}/data"
    rm "${SCRIPT_DIR}/data/micrococcus_data.tar.gz"
fi

# Run R analysis (staying in the same environment since R is already installed)
echo "Running comparative genomics analysis..."
Rscript run_comparative.R

echo "Analysis complete! Check the results in the outputs/ and results/ directories."