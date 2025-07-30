#!/bin/bash
set -e

# Set dynamic thread count
THREADS=$(nproc)

# Logging functions
TOTAL_STEPS=8
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

# STEP 1: Environment Setup
log_step "Setting up conda environments"

if ! env_exists "metagenomics"; then
    log_substep "Creating metagenomics environment"
    mamba env create -f environment.yml -y
fi

log_completion

# Create directory structure
mkdir -p data
mkdir -p outputs/0_fastqc
mkdir -p outputs/1_trimmed
mkdir -p outputs/2_fastqc
mkdir -p outputs/3_assembly
mkdir -p outputs/4_taxonomy
mkdir -p outputs/5_biom
mkdir -p results
mkdir -p data/reference

# STEP 2: Data Download
log_step "Downloading metagenomics data"

log_substep "Downloading data"
if [ ! -f "data/metagenomics.zip" ]; then
    wget -O data/metagenomics.zip https://zenodo.org/records/7010950/files/dc_workshop.zip
fi

log_substep "Downloading metadata"
if [ ! -f "data/MGRAST_MetaData_JP.xlsx" ]; then
    wget -O data/MGRAST_MetaData_JP.xlsx https://zenodo.org/records/7010950/files/MGRAST_MetaData_JP.xlsx
fi

log_substep "Downloading Kraken2 standard database"
if [ ! -f "data/k2_standard_16gb_20241228.tar.gz" ]; then
    wget -v -O data/k2_standard_16gb_20241228.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz
fi

log_substep "Extracting data"
unzip data/metagenomics.zip -d data/

log_substep "Moving extracted files to expected locations"
mv data/dc_workshop/data/untrimmed_fastq data/

log_substep "Extracting Kraken2 database"
if [ ! -d "data/reference/hash.k2d" ]; then
    tar -xvf data/k2_standard_16gb_20241228.tar.gz -C data/reference
fi

log_completion

# STEP 3: Initial Quality Control
log_step "Initial quality control with FastQC"

log_substep "Running FastQC on raw reads"
mamba run -n metagenomics fastqc data/untrimmed_fastq/*fastq.gz -o outputs/0_fastqc -t $THREADS

log_completion

# STEP 4: Read Trimming
log_step "Trimming reads with Trimmomatic"

log_substep "Decompressing FASTQ files"
for file in data/untrimmed_fastq/*.fastq.gz; do
    gunzip -k $file
done

log_substep "Trimming paired-end reads"
for file in data/untrimmed_fastq/*_R1.fastq; do
    base=$(basename "$file" _R1.fastq) 
    log_substep "Processing sample: $base"
    mamba run -n metagenomics trimmomatic PE -threads $THREADS \
        data/untrimmed_fastq/${base}_R1.fastq \
        data/untrimmed_fastq/${base}_R2.fastq \
        outputs/1_trimmed/${base}_1P.fastq \
        outputs/1_trimmed/${base}_1U.fastq \
        outputs/1_trimmed/${base}_2P.fastq \
        outputs/1_trimmed/${base}_2U.fastq \
        SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:data/untrimmed_fastq/TruSeq3-PE.fa:2:40:15
done

log_completion

# STEP 5: Post-trimming Quality Control
log_step "Post-trimming quality control"

log_substep "Running FastQC on trimmed reads"
mamba run -n metagenomics fastqc outputs/1_trimmed/*.fastq -o outputs/2_fastqc -t $THREADS

log_substep "Generating MultiQC report"
mamba run -n metagenomics multiqc outputs/2_fastqc -o outputs/2_fastqc

log_completion

# STEP 6: Metagenomic Assembly
log_step "Metagenomic assembly with SPAdes"

log_substep "Assembling JC1A sample"
mamba run -n metagenomics metaspades.py -1 outputs/1_trimmed/JC1A_1P.fastq \
                                      -2 outputs/1_trimmed/JC1A_2P.fastq \
                                      -o outputs/3_assembly/assembly_JC1A \
                                      -t $THREADS

log_substep "Assembling JP4D sample"
mamba run -n metagenomics metaspades.py -1 outputs/1_trimmed/JP4D_1P.fastq \
                                      -2 outputs/1_trimmed/JP4D_2P.fastq \
                                      -o outputs/3_assembly/assembly_JP4D \
                                      -t $THREADS

log_completion

# STEP 7: Taxonomic Classification
log_step "Taxonomic classification with Kraken2"

log_substep "Classifying JC1A assembly"
mamba run -n metagenomics kraken2 --db data/reference --threads $THREADS \
                               outputs/3_assembly/assembly_JC1A/contigs.fasta \
                               --output "outputs/4_taxonomy/JC1A.kraken" \
                               --report "outputs/4_taxonomy/JC1A.report"

log_substep "Classifying JP4D assembly"
mamba run -n metagenomics kraken2 --db data/reference --threads $THREADS \
                               outputs/3_assembly/assembly_JP4D/contigs.fasta \
                               --output "outputs/4_taxonomy/JP4D.kraken" \
                               --report "outputs/4_taxonomy/JP4D.report"

log_substep "Converting to BIOM format"
mamba run -n metagenomics kraken-biom outputs/4_taxonomy/JC1A.report \
                                   outputs/4_taxonomy/JP4D.report \
                                   --fmt json -o outputs/5_biom/cuatroc.biom

log_completion

# STEP 8: Phyloseq Analysis
log_step "Running phyloseq analysis"

log_substep "Executing R phyloseq script"
mamba run -n metagenomics Rscript run_phyloseq.R

log_completion

echo "=================================================================="
echo "METAGENOMICS PIPELINE COMPLETED SUCCESSFULLY"
echo "=================================================================="
echo "Results are available in results/"
echo "Completed at: $(date)"