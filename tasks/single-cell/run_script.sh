#!/bin/bash
set -e

mkdir -p ./data
mkdir -p ./reference
mkdir -p ./results

echo "Downloading scRNA-seq data files..."

declare -a samples=(
    "GSM6611295:GSM6611295_P15306_5001"
    "GSM6611296:GSM6611296_P15306_5002"
    "GSM6611297:GSM6611297_P14601_4004"
    "GSM6611298:GSM6611298_P14601_4005"
    "GSM6611299:GSM6611299_P15306_5003"
    "GSM6611300:GSM6611300_P15306_5004"
)

download_and_decompress() {
    local accession=$1
    local file_prefix=$2
    local file_type=$3
    local ext=$4
    
    local gz_file="./data/${file_prefix}_${file_type}.${ext}.gz"
    local target_file="./data/${file_prefix}_${file_type}.${ext}"
    
    if [ -f "$target_file" ]; then
        echo "File $target_file already exists, skipping download"
        return
    fi
    
    if [ ! -f "$gz_file" ]; then
        local encoded_prefix=$(echo "$file_prefix" | sed 's/_/%5F/g' | sed 's/\./%2E/g')
        local url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${accession}&format=file&file=${encoded_prefix}%5F${file_type}.${ext}.gz"
        
        echo "Downloading ${file_type} file for ${accession}..."
        curl -L -o "$gz_file" "$url"
        
        if [ $? -ne 0 ]; then
            echo "Error downloading $gz_file"
            exit 1
        fi
    fi
    
    echo "Decompressing $gz_file..."
    gunzip -f "$gz_file"
}

for sample in "${samples[@]}"; do
    IFS=':' read -r accession file_prefix <<< "$sample"
    
    echo "Processing sample: $accession ($file_prefix)"
    
    download_and_decompress "$accession" "$file_prefix" "matrix" "mtx"
    download_and_decompress "$accession" "$file_prefix" "barcodes" "tsv"
    download_and_decompress "$accession" "$file_prefix" "features" "tsv"
done

echo "Data download completed!"

echo "Downloading cell marker database..."
if [ ! -f "./reference/Cell_marker_Seq.xlsx" ]; then
    curl -L -o "./reference/Cell_marker_Seq.xlsx" "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Seq.xlsx"
    if [ $? -ne 0 ]; then
        echo "Error downloading cell marker database"
        exit 1
    fi
    echo "Cell marker database downloaded successfully!"
else
    echo "Cell marker database already exists, skipping download"
fi

python run_analysis.py