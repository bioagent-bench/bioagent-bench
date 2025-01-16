#!/bin/bash
set -e

input_dir="/home/dev/benchmark/bio-agent-benchmark/scan/data/cells-dataset/fast-dump"
output_dir="/home/dev/benchmark/bio-agent-benchmark/scan/data/cells-dataset/filtered"
report_dir="/home/dev/benchmark/bio-agent-benchmark/scan/data/cells-dataset/qc_reports"

mkdir -p "$output_dir" "$report_dir"

process_fastq_pair() {
    local input_file="$1"
    local base_name=$(basename "$input_file" _1.fastq)
    echo "Processing: $base_name"
    
    /home/dev/benchmark/bio-agent-benchmark/scan/scripts/fastp \
        -i "$input_file" \
        -I "$(dirname "$input_file")/${base_name}_2.fastq" \
        -o "$output_dir/${base_name}_1.filtered.fastq.gz" \
        -O "$output_dir/${base_name}_2.filtered.fastq.gz" \
        --json "$report_dir/${base_name}_fastp.json" \
        --thread 1 \
        2> "$report_dir/${base_name}_fastp.log" || {
            echo "Error processing $base_name" >&2
            return 1
        }
}
export -f process_fastq_pair
export input_dir output_dir report_dir

if [ ! -d "$input_dir" ]; then
    echo "Input directory does not exist: $input_dir"
    exit 1
fi

# find all *_1.fastq files and process them in parallel
find "$input_dir" -name "*_1.fastq" | \
    parallel --jobs $(nproc) process_fastq_pair {}

echo "All quality control processing completed"
