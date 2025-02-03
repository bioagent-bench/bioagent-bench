#!/bin/bash

DATA_PATH="/home/dev/benchmark/bio-agent-benchmark/virus/data"
OUTPUT_PATH="/home/dev/benchmark/bio-agent-benchmark/virus/data/fastq_output"
mkdir -p "$OUTPUT_PATH"

for sra_file in "$DATA_PATH"/SRR*/*.sra; do
    if [ -e "$sra_file" ]; then
        echo "Processing $sra_file..."
        fasterq-dump "$sra_file" --split-files -O "$OUTPUT_PATH"
    else
        echo "No SRA files found in $DATA_PATH/SRR*"
        break
    fi
done

echo "All SRA files processed."
