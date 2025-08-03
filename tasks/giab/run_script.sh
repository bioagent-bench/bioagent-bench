#!/bin/bash

set -euo pipefail 

# Set dynamic thread count
THREADS=$(nproc)

# Logging functions
TOTAL_STEPS=13
CURRENT_STEP=0

log_step() {
    CURRENT_STEP=$((CURRENT_STEP + 1))
    echo "=================================================================="
    echo "STEP $CURRENT_STEP/$TOTAL_STEPS: $1"
    echo "=================================================================="
    echo "Started at: $(date)"
    echo ""
}

log_substep() {
    echo "  → $1"
}

log_completion() {
    echo ""
    echo "Completed at: $(date)"
    echo ""
}

env_exists() {
    local env_name=$1
    mamba env list | grep -q "^${env_name} "
}

echo "Creating data directory..."
mkdir -p ./benchmark_data
mkdir -p ./data
mkdir -p ./fastq_data
mkdir -p ./benchmark_ref/
mkdir -p ./fastq_merged/
mkdir -p ./bam/

BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38"
FILES=(
    "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
    "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
)

echo "Starting download of GIAB benchmark files..."

# Download each file
for file in "${FILES[@]}"; do
    echo "Downloading $file..."

    if [[ -f "./benchmark_data/$file" ]]; then
        echo "  $file already exists, skipping..."
        continue
    fi
    wget -O "./benchmark_data/$file" "$BASE_URL/$file"
done


echo ""
echo "Download complete! Files saved to ./benchmark_data/"
echo "Downloaded files:"
ls -lh ./benchmark_data/

if ! env_exists "giab"; then
    log_substep "Creating giab environment"
    mamba env create -f environment.yml
fi

echo ""
echo "Starting chromosome 20 filtering..."

# Subset VCF to chromosome 20
echo "Subsetting VCF to chr20..."
if [[ -f "./benchmark_data/HG001_chr20_truth.vcf.gz" ]]; then
    echo "  HG001_chr20_truth.vcf.gz already exists, skipping..."
else
    mamba run -n giab bcftools view -r chr20 ./benchmark_data/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -Oz -o ./benchmark_data/HG001_chr20_truth.vcf.gz
    # Index the output VCF
    mamba run -n giab tabix -p vcf ./benchmark_data/HG001_chr20_truth.vcf.gz
    echo "  Created HG001_chr20_truth.vcf.gz"
fi

# Subset BED to chromosome 20
echo "Subsetting BED to chr20..."
if [[ -f "./benchmark_data/HG001_chr20_highconf.bed" ]]; then
    echo "  HG001_chr20_highconf.bed already exists, skipping..."
else
    grep -w "^chr20" ./benchmark_data/HG001_GRCh38_1_22_v4.2.1_benchmark.bed > ./benchmark_data/HG001_chr20_highconf.bed
    echo "  Created HG001_chr20_highconf.bed"
fi

echo ""
echo "Filtering complete! Chromosome 20 files:"
echo "Truth VCF: ./benchmark_data/HG001_chr20_truth.vcf.gz"
echo "High-confidence regions BED: ./benchmark_data/HG001_chr20_highconf.bed"
echo ""
echo "File details:"
ls -lh ./benchmark_data/HG001_chr20_*

echo ""
echo "Starting download of FASTQ files..."

FASTQ_BASE_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0a"
FASTQ_FILES=(
    "SampleSheet.csv"
    "U0a_CGATGT_L001_R1_001.fastq.gz"
    "U0a_CGATGT_L001_R1_002.fastq.gz"
    "U0a_CGATGT_L001_R1_003.fastq.gz"
    "U0a_CGATGT_L001_R1_004.fastq.gz"
    "U0a_CGATGT_L001_R1_005.fastq.gz"
    "U0a_CGATGT_L001_R2_001.fastq.gz"
    "U0a_CGATGT_L001_R2_002.fastq.gz"
    "U0a_CGATGT_L001_R2_003.fastq.gz"
    "U0a_CGATGT_L001_R2_004.fastq.gz"
    "U0a_CGATGT_L001_R2_005.fastq.gz"
    "U0a_CGATGT_L002_R1_001.fastq.gz"
    "U0a_CGATGT_L002_R1_002.fastq.gz"
    "U0a_CGATGT_L002_R1_003.fastq.gz"
    "U0a_CGATGT_L002_R1_004.fastq.gz"
    "U0a_CGATGT_L002_R1_005.fastq.gz"
    "U0a_CGATGT_L002_R2_001.fastq.gz"
    "U0a_CGATGT_L002_R2_002.fastq.gz"
    "U0a_CGATGT_L002_R2_003.fastq.gz"
    "U0a_CGATGT_L002_R2_004.fastq.gz"
    "U0a_CGATGT_L002_R2_005.fastq.gz"
)

all_files_exist=true
for file in "${FASTQ_FILES[@]}"; do
    if [[ ! -f "./fastq_data/$file" ]]; then
        all_files_exist=false
        break
    fi
done

if [[ "$all_files_exist" == "true" ]]; then
    echo "All FASTQ files already exist in ./fastq_data/, skipping download..."
else
    echo "Downloading all FASTQ files from GIAB repository..."
    echo "This may take a while"
    echo "Downloading ${#FASTQ_FILES[@]} files..."
    
    # Download each file individually
    for i in "${!FASTQ_FILES[@]}"; do
        file="${FASTQ_FILES[$i]}"
        echo "Downloading file $((i+1))/${#FASTQ_FILES[@]}: $file"
        
        if [[ -f "./fastq_data/$file" ]]; then
            echo "  $file already exists, skipping..."
            continue
        fi
        
        wget -O "./fastq_data/$file" "$FASTQ_BASE_URL/$file"
    done
    
    echo "FASTQ download complete!"
fi

echo ""
echo "FASTQ files downloaded to ./fastq_data/"
echo "Downloaded files:"
ls -lh ./fastq_data/

echo ""
echo "Starting download of reference files..."

wget -O ./benchmark_ref/GRCh38_noalt.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz

echo "Reference file downloaded to ./benchmark_ref/GRCh38_noalt.fa.gz"

echo "Unzipping reference file..."
gunzip ./benchmark_ref/GRCh38_noalt.fa.gz

echo "Reference file unzipped to ./benchmark_ref/GRCh38_noalt.fa"

echo "Indexing reference file..."
mamba run -n giab bwa index -a bwtsw -p ./benchmark_ref/GRCh38_noalt ./benchmark_ref/GRCh38_noalt.fa
echo "Indexing complete!"

echo "Merging FASTQ files..."
cat fastq_data/U0a_CGATGT_L00{1..2}_R1_*.fastq.gz > fastq_merged/U0a_R1.fq.gz
cat fastq_data/U0a_CGATGT_L00{1..2}_R2_*.fastq.gz > fastq_merged/U0a_R2.fq.gz
echo "Merging complete!"

bwa mem -t 64 \
        -R '@RG\tID:U0a\tSM:NA12878\tPL:ILLUMINA' \
        ./benchmark_ref/GRCh38_noalt \
        fastq_merged/U0a_R1.fq.gz fastq_merged/U0a_R2.fq.gz  \
    | samtools sort -@64 -o ./bam/HG001_U0a.all.bam

samtools index bam/HG001_U0a.all.bam

echo "Subsetting BAM file to chr20..."
samtools view -@64 -b -h bam/HG001_U0a.all.bam chr20 > bam/HG001_U0a.chr20.bam
samtools index bam/HG001_U0a.chr20.bam
samtools idxstats bam/HG001_U0a.chr20.bam

echo "Converting back to FastQ..."
samtools fastq -@64                         \
    -1 data/HG001_chr20_R1.fq.gz            \
    -2 data/HG001_chr20_R2.fq.gz            \
    -0 /dev/null -s /dev/null -n           \
    bam/HG001_U0a.chr20.bam
echo "Conversion complete!"






