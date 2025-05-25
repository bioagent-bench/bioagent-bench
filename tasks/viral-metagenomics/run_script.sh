#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Calculate optimal thread count (75% of available cores)
THREADS=$(( $(nproc) * 3 / 4 ))

source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

mamba env create -f "${SCRIPT_DIR}/environment.yml" || true
mamba activate viral-metagenomics

# Create directory structure
mkdir -p "${SCRIPT_DIR}/data"
mkdir -p "${SCRIPT_DIR}/outputs"

# Download input data
wget -O "${SCRIPT_DIR}/data/Dol1_S19_L001_R1_001.fastq.gz" https://osf.io/4x6qs/download
wget -O "${SCRIPT_DIR}/data/Dol1_S19_L001_R2_001.fastq.gz" https://osf.io/z2xed/download

# Download reference data
mkdir -p "${SCRIPT_DIR}/data/reference"
wget -O "${SCRIPT_DIR}/data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz" https://ftp.ensembl.org/pub/release-113/fasta/tursiops_truncatus/dna/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz

# 1. Trimming
mkdir -p "${SCRIPT_DIR}/outputs/1_trimming"
fastp -i "${SCRIPT_DIR}/data/Dol1_S19_L001_R1_001.fastq.gz" \
      -o "${SCRIPT_DIR}/outputs/1_trimming/Dol1_trimmed_R1.fastq" \
      -I "${SCRIPT_DIR}/data/Dol1_S19_L001_R2_001.fastq.gz" \
      -O "${SCRIPT_DIR}/outputs/1_trimming/Dol1_trimmed_R2.fastq" \
      --detect_adapter_for_pe --length_required 30 \
      --cut_front --cut_tail --cut_mean_quality 10 \
      --thread ${THREADS}

# Prepare reference
gunzip -k "${SCRIPT_DIR}/data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz"

# 2. Indexing and 3. Mapping
mkdir -p "${SCRIPT_DIR}/outputs/2_indexing"
mkdir -p "${SCRIPT_DIR}/outputs/3_mapping"

bowtie2-build "${SCRIPT_DIR}/data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa" \
              "${SCRIPT_DIR}/outputs/2_indexing/Tursiops_truncatus" \
              --threads ${THREADS}

bowtie2 -x "${SCRIPT_DIR}/outputs/2_indexing/Tursiops_truncatus" \
        -1 "${SCRIPT_DIR}/outputs/1_trimming/Dol1_trimmed_R1.fastq" \
        -2 "${SCRIPT_DIR}/outputs/1_trimming/Dol1_trimmed_R2.fastq" \
        -S "${SCRIPT_DIR}/outputs/3_mapping/dol_map.sam" \
        --un-conc "${SCRIPT_DIR}/outputs/3_mapping/Dol_reads_unmapped.fastq" \
        --threads ${THREADS}

# 4. Assembly
mkdir -p "${SCRIPT_DIR}/outputs/4_assembly"
megahit -1 "${SCRIPT_DIR}/outputs/3_mapping/Dol_reads_unmapped.1.fastq" \
        -2 "${SCRIPT_DIR}/outputs/3_mapping/Dol_reads_unmapped.2.fastq" \
        --num-cpu-threads ${THREADS} \
        -o "${SCRIPT_DIR}/outputs/4_assembly/megahit_out"

# 5. Classification
mkdir -p "${SCRIPT_DIR}/outputs/5_classification"
cd "${SCRIPT_DIR}/outputs/5_classification"
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-08-15.tgz
tar xzf kaiju_db_viruses_2024-08-15.tgz
rm kaiju_db_viruses_2024-08-15.tgz

kaiju -t "${SCRIPT_DIR}/outputs/5_classification/nodes.dmp" \
      -f "${SCRIPT_DIR}/outputs/5_classification/kaiju_db_viruses.fmi" \
      -i "${SCRIPT_DIR}/outputs/4_assembly/megahit_out/final.contigs.fa" \
      -o "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.out" \
      -z ${THREADS}

kaiju2krona -t "${SCRIPT_DIR}/outputs/5_classification/nodes.dmp" \
            -n "${SCRIPT_DIR}/outputs/5_classification/names.dmp" \
            -i "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.out" \
            -o "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.krona" \
            -u

ktImportText -o "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.krona.html" \
             "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.krona"

# Results
mkdir -p "${SCRIPT_DIR}/results"

# Create viral taxonomy results
echo "contig_count,domain,species" > "${SCRIPT_DIR}/results/taxonomy.csv"
awk -F'\t' 'NF>1 {
    species = $(NF);  # Last column
    if (NF == 1) species = "NA";  # Handle Unclassified case
    print $1 "," $3 "," species
}' "${SCRIPT_DIR}/outputs/5_classification/Dol1_contigs_kaiju.krona" >> "${SCRIPT_DIR}/results/taxonomy.csv"