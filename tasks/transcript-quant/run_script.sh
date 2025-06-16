#!/bin/bash

# Exit on error
set -e

# Get number of available CPU threads
THREADS=$(nproc)

# Create necessary directories
mkdir -p processing/index
mkdir -p processing/quant
mkdir -p results

# Define input/output paths
TRANSCRIPTS="data/transcriptome.fa"
INDEX_DIR="processing/index"
QUANT_DIR="processing/quant"
READS1="data/reads_1.fq.gz"
READS2="data/reads_2.fq.gz"
FINAL_RESULTS="results/results.tsv"

# Run Salmon index
echo "Running Salmon index..."
salmon index -t "$TRANSCRIPTS" -i "$INDEX_DIR" 2>&1 | tee "$INDEX_DIR/index.log"

# Run Salmon quant
echo "Running Salmon quant..."
salmon quant \
    -i "$INDEX_DIR" \
    -l A \
    -p "$THREADS" \
    -1 "$READS1" \
    -2 "$READS2" \
    --validateMappings \
    -o "$QUANT_DIR" 2>&1 | tee "$QUANT_DIR/quant.log"

# Process the quant.sf file to match the truth.tsv format
echo "Processing results..."
awk 'NR>1 {print $1 "\t" int($5)}' "$QUANT_DIR/quant.sf" > "$FINAL_RESULTS"

echo "Analysis complete! Results saved to $FINAL_RESULTS"