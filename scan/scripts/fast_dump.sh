#!/bin/bash
base_dir="/home/dev/benchmark/scan/data/cells-dataset/SRR9136455"

# Check if base directory exists
if [ ! -d "$base_dir" ]; then
    echo "Base directory does not exist: $base_dir"
    exit 1
fi

# Function to process a single SRA file
process_sra() {
    local sra_file="$1"
    local dir="$2"
    echo "Converting file: $sra_file in $dir"
    cd "$dir" || return 1
    
    if fastq-dump --split-files "$sra_file" -O /home/dev/benchmark/scan/data/cells_dataset/fast-dump; then
        echo "Successfully converted: $sra_file"
    else
        echo "Error converting: $sra_file"
        return 1
    fi
}
export -f process_sra

# Find all SRA files and process them in parallel
find "$base_dir" -type f -name "*.sra" | \
    parallel --jobs $(nproc) process_sra {} {//}

echo "All conversions completed"
