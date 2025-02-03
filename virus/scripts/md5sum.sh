#!/bin/bash

# Relative path from the scripts directory to the data/fastq_output folder.
FASTQ_DIR="../data/fastq_output"
OUTPUT_MD5="../data/fastq_output/md5sum.txt"
SAMPLESHEET="../data/fastq_output/samplesheet.tsv"

# Remove any previous output files.
rm -f "$OUTPUT_MD5" "$SAMPLESHEET"

# Write the header for the samplesheet TSV file.
echo -e "sample_id\tpair_id\tfraction\tR1\tR2\tR1_MD5\tR2_MD5" > "$SAMPLESHEET"

# Predefined sample data from the original paper.
read -r -d '' SAMPLE_DATA << 'EOF'
SRR8487026	R1	microbial
SRR8487027	R2	microbial
SRR8487030	R3	microbial
SRR8487031	R4	microbial
SRR8487039	R5	microbial
SRR8487040	R6	microbial
SRR8487012	R7	microbial
SRR8487013	R8	microbial
SRR8487035	R1	viral
SRR8487036	R2	viral
SRR8487015	R3	viral
SRR8487014	R4	viral
SRR8487018	R5	viral
SRR8487021	R6	viral
SRR8487034	R7	viral
SRR8487023	R8	viral
EOF

# Process each sample line from the predefined data.
while IFS=$'\t' read -r sample_id pair_id fraction; do
    echo "--------------------------------------------------"
    echo "Processing sample: ${sample_id} (${fraction}, pair_id: ${pair_id})"

    # Construct paths to the FASTQ files (uncompressed)
    r1_file="${FASTQ_DIR}/${sample_id}_1.fastq"
    r2_file="${FASTQ_DIR}/${sample_id}_2.fastq"

    if [ ! -f "$r1_file" ]; then
        echo "Warning: $r1_file not found."
    fi
    if [ ! -f "$r2_file" ]; then
        echo "Warning: $r2_file not found."
    fi

    # Calculate MD5 checksums.
    r1_md5=$(md5sum "$r1_file" 2>/dev/null | awk '{print $1}')
    r2_md5=$(md5sum "$r2_file" 2>/dev/null | awk '{print $1}')

    echo "Checksum for R1 ($r1_file): $r1_md5"
    echo "Checksum for R2 ($r2_file): $r2_md5"

    # Append the MD5 information to the MD5 output file.
    echo "$r1_md5  $r1_file" >> "$OUTPUT_MD5"
    echo "$r2_md5  $r2_file" >> "$OUTPUT_MD5"

    # Append the sample information along with the MD5 checksums to the samplesheet.
    echo -e "${sample_id}\t${pair_id}\t${fraction}\t${r1_file}\t${r2_file}\t${r1_md5}\t${r2_md5}" >> "$SAMPLESHEET"

    echo "Finished processing sample: ${sample_id}"
done <<< "$SAMPLE_DATA"

echo "--------------------------------------------------"
echo "MD5 checksums have been saved to: $OUTPUT_MD5"
echo "Samplesheet TSV has been generated at: $SAMPLESHEET"
