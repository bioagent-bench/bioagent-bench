#!/bin/bash

FASTQ_DIR=/home/dev/benchmark/bio-agent-benchmark/virus/data/fastq_output
OUTPUT_FILE=/home/dev/benchmark/bio-agent-benchmark/virus/data/fastq_output/md5sum.txt

if [ -f "$OUTPUT_FILE" ]; then
    rm "$OUTPUT_FILE"
fi

# Loop through all the FASTQ files and calculate MD5 for each.
for file in "$FASTQ_DIR"/*.fastq; do
    echo "Processing $file..."
    md5sum "$file" >> "$OUTPUT_FILE"
done

echo "MD5 checksums have been saved to $OUTPUT_FILE"
