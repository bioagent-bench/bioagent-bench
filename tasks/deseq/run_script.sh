#!/bin/bash
set -e

# Get the absolute path of the script directory
SCRIPT_DIR="/app"
# Calculate optimal thread count (75% of available cores)
THREADS=$(( $(nproc) * 3 / 4 ))

# Activate the main environment
# source ~/.bashrc
# mamba activate rnaseq

# Create output directories
mkdir -p "${SCRIPT_DIR}/data/processing"/{0_fasterqdump,1_fastqc,2_multiqc,3_trimming,4_fastqc_trimmed,5_multiqc_trimmed,6_indexing,7_mapping,8_deseq,9_gorich}
mkdir -p "${SCRIPT_DIR}/data/reference"

# Download SRA files if they don't exist
prefetch SRR1278968 SRR1278969 SRR1278970 SRR1278971 SRR1278972 SRR1278973 \
    -O "${SCRIPT_DIR}/data/"

# Convert SRA to FASTQ
echo "Converting SRA to FASTQ..."
for sra_file in "${SCRIPT_DIR}"/data/SRR*/*.sra; do
    fasterq-dump "$sra_file" -O "${SCRIPT_DIR}/data/processing/0_fasterqdump"
done 

# Run FastQC on raw data
echo "Running FastQC on raw data..."
fastqc "${SCRIPT_DIR}/data/processing/0_fasterqdump/"*fastq \
    -o "${SCRIPT_DIR}/data/processing/1_fastqc" \
    -t "${THREADS}"
multiqc "${SCRIPT_DIR}/data/processing/1_fastqc" \
    -o "${SCRIPT_DIR}/data/processing/2_multiqc"

# Trim reads
echo "Trimming reads..."
for file in "${SCRIPT_DIR}/data/processing/0_fasterqdump/"*_1.fastq; do
    base=$(basename "$file" _1.fastq) 
    trimmomatic PE -threads "${THREADS}" \
        "${SCRIPT_DIR}/data/processing/0_fasterqdump/${base}_1.fastq" \
        "${SCRIPT_DIR}/data/processing/0_fasterqdump/${base}_2.fastq" \
        "${SCRIPT_DIR}/data/processing/3_trimming/${base}_1P.fastq" \
        "${SCRIPT_DIR}/data/processing/3_trimming/${base}_1U.fastq" \
        "${SCRIPT_DIR}/data/processing/3_trimming/${base}_2P.fastq" \
        "${SCRIPT_DIR}/data/processing/3_trimming/${base}_2U.fastq" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done

# Run FastQC on trimmed data
echo "Running FastQC on trimmed data..."
fastqc "${SCRIPT_DIR}/data/processing/3_trimming/"*fastq \
    -o "${SCRIPT_DIR}/data/processing/4_fastqc_trimmed" \
    -t "${THREADS}"
multiqc "${SCRIPT_DIR}/data/processing/4_fastqc_trimmed" \
    -o "${SCRIPT_DIR}/data/processing/5_multiqc_trimmed"

# Download reference data if it doesn't exist
echo "Downloading reference data..."
if [ ! -f "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta" ]; then
    wget -O "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz" \
        http://www.candidagenome.org/download/sequence/C_parapsilosis_CDC317/current/C_parapsilosis_CDC317_current_chromosomes.fasta.gz
    gunzip "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz"
fi

if [ ! -f "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_features.gff" ]; then
    wget -O "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_features.gff" \
        http://www.candidagenome.org/download/gff/C_parapsilosis_CDC317/C_parapsilosis_CDC317_current_features.gff
fi

# Generate STAR index
echo "Generating STAR index..."
STAR --runThreadN "${THREADS}" \
     --runMode genomeGenerate \
     --genomeDir "${SCRIPT_DIR}/data/processing/6_indexing/" \
     --genomeFastaFiles "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta" \
     --genomeSAindexNbases 10

# Convert GFF to GTF
echo "Converting GFF to GTF..."
gffread "${SCRIPT_DIR}/data/reference/C_parapsilosis_CDC317_current_features.gff" \
    -T -o "${SCRIPT_DIR}/data/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf"

# Run STAR mapping
echo "Running STAR mapping..."
for file in "${SCRIPT_DIR}/data/processing/0_fasterqdump/"*_1.fastq; do
    base=$(basename "$file" _1.fastq)
    STAR --runThreadN "${THREADS}" \
        --genomeDir "${SCRIPT_DIR}/data/processing/6_indexing" \
        --sjdbGTFfile "${SCRIPT_DIR}/data/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf" \
        --readFilesIn "${SCRIPT_DIR}/data/processing/3_trimming/${base}_1P.fastq" \
                      "${SCRIPT_DIR}/data/processing/3_trimming/${base}_2P.fastq" \
        --outFileNamePrefix "${SCRIPT_DIR}/data/processing/7_mapping/${base}_" \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 100000000000 \
        --quantMode GeneCounts
done

# Generate mapping QC report
echo "Generating mapping QC report..."
multiqc "${SCRIPT_DIR}/data/processing/7_mapping" \
    -o "${SCRIPT_DIR}/data/processing/7_mapping"

# Download GO annotations
echo "Downloading GO annotations..."
if [ ! -f "${SCRIPT_DIR}/data/processing/9_gorich/gene_association.cgd" ]; then
    wget -O "${SCRIPT_DIR}/data/processing/9_gorich/gene_association.cgd.gz" \
        http://www.candidagenome.org/download/go/gene_association.cgd.gz
    gunzip "${SCRIPT_DIR}/data/processing/9_gorich/gene_association.cgd.gz"
fi

# Run R analysis
mkdir -p results
python run_deseq.py

echo "Analysis complete! Check the results in the data/processing directory."