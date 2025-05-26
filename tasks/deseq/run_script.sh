#!/bin/bash
set -e

# Get the absolute path of the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Calculate optimal thread count (75% of available cores)
THREADS=$(( $(nproc) * 3 / 4 ))

# Source conda/mamba
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# Create and activate the main environment
mamba env create -f "${SCRIPT_DIR}/environment.yml" || true
mamba activate rnaseq

# Create output directories
mkdir -p "${SCRIPT_DIR}/processing"/{0_fasterqdump,1_fastqc,2_multiqc,3_trimming,4_fastqc_trimmed,5_multiqc_trimmed,6_indexing,7_mapping,8_deseq,9_gorich}
mkdir -p "${SCRIPT_DIR}/data/reference"

# Download SRA files
# prefetch SRR1278968 SRR1278969 SRR1278970 SRR1278971 SRR1278972 SRR1278973 \
#     -O "${SCRIPT_DIR}/data/"

# Convert SRA to FASTQ
for sra_file in "${SCRIPT_DIR}"/data/SRR*/*.sra; do
    fasterq-dump "$sra_file" -O "${SCRIPT_DIR}/processing/0_fasterqdump"
done 

# Run FastQC on raw data
fastqc "${SCRIPT_DIR}/processing/0_fasterqdump/"*fastq \
    -o "${SCRIPT_DIR}/processing/1_fastqc" \
    -t "${THREADS}"
multiqc "${SCRIPT_DIR}/processing/1_fastqc" \
    -o "${SCRIPT_DIR}/processing/2_multiqc"

# Trim reads
for file in "${SCRIPT_DIR}/processing/0_fasterqdump/"*_1.fastq; do
    base=$(basename "$file" _1.fastq) 
    trimmomatic PE -threads "${THREADS}" \
        "${SCRIPT_DIR}/processing/0_fasterqdump/${base}_1.fastq" \
        "${SCRIPT_DIR}/processing/0_fasterqdump/${base}_2.fastq" \
        "${SCRIPT_DIR}/processing/3_trimming/${base}_1P.fastq" \
        "${SCRIPT_DIR}/processing/3_trimming/${base}_1U.fastq" \
        "${SCRIPT_DIR}/processing/3_trimming/${base}_2P.fastq" \
        "${SCRIPT_DIR}/processing/3_trimming/${base}_2U.fastq" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done

# Run FastQC on trimmed data
fastqc "${SCRIPT_DIR}/processing/3_trimming/"*fastq \
    -o "${SCRIPT_DIR}/processing/4_fastqc_trimmed" \
    -t "${THREADS}"
multiqc "${SCRIPT_DIR}/processing/4_fastqc_trimmed" \
    -o "${SCRIPT_DIR}/processing/5_multiqc_trimmed"

# Download reference data
wget -O "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz" \
    http://www.candidagenome.org/download/sequence/C_parapsilosis_CDC317/current/C_parapsilosis_CDC317_current_chromosomes.fasta.gz
wget -O "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_features.gff" \
    http://www.candidagenome.org/download/gff/C_parapsilosis_CDC317/C_parapsilosis_CDC317_current_features.gff
gunzip -k "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz"

# Generate STAR index
STAR --runThreadN "${THREADS}" \
     --runMode genomeGenerate \
     --genomeDir "${SCRIPT_DIR}/processing/6_indexing/" \
     --genomeFastaFiles "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta" \
     --genomeSAindexNbases 10

# Convert GFF to GTF
gffread "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_features.gff" \
    -T -o "${SCRIPT_DIR}/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf"

# Run STAR mapping
for file in "${SCRIPT_DIR}/processing/0_fasterqdump/"*_1.fastq; do
    base=$(basename "$file" _1.fastq)
    STAR --runThreadN "${THREADS}" \
        --genomeDir "${SCRIPT_DIR}/processing/6_indexing" \
        --sjdbGTFfile "${SCRIPT_DIR}/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf" \
        --readFilesIn "${SCRIPT_DIR}/processing/3_trimming/${base}_1P.fastq" \
                      "${SCRIPT_DIR}/processing/3_trimming/${base}_2P.fastq" \
        --outFileNamePrefix "${SCRIPT_DIR}/processing/7_mapping/${base}_" \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 100000000000 \
        --quantMode GeneCounts
done

# Generate mapping QC report
multiqc "${SCRIPT_DIR}/processing/7_mapping" \
    -o "${SCRIPT_DIR}/processing/7_mapping"

# Download GO annotations
wget -O "${SCRIPT_DIR}/processing/9_gorich/gene_association.cgd.gz" \
    http://www.candidagenome.org/download/go/gene_association.cgd.gz
gunzip -k "${SCRIPT_DIR}/processing/9_gorich/gene_association.cgd.gz"

# Switch to R environment and run analysis
mamba deactivate
mamba env create -f "${SCRIPT_DIR}/r-environment.yml" || true
mamba activate rnaseq-r

# Run R analysis
Rscript "${SCRIPT_DIR}/run_deseq.R"