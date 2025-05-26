#!/bin/bash
set -e

# Get the absolute path of the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Calculate optimal memory (70% of total system memory)
TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
AVAILABLE_MEM_GB=$((TOTAL_MEM_KB * 7 / (1024 * 1024 * 10)))

source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# Create directory structure
mkdir -p "${SCRIPT_DIR}/data"

# Download and extract source data
wget -O "${SCRIPT_DIR}/data/protocols.zip" \
    "http://sourceforge.net/projects/snpeff/files/protocols.zip"
unzip "${SCRIPT_DIR}/data/protocols.zip" -d "${SCRIPT_DIR}/data/"
rm "${SCRIPT_DIR}/data/protocols.zip"  # Clean up zip file after extraction

# Set up environment
mamba env create -f "${SCRIPT_DIR}/environment.yml" || true
mamba activate fibrosis_env

# Download reference data
snpEff download -v GRCh37.75

# Run SNP effect prediction
snpEff \
    -Xmx${AVAILABLE_MEM_GB}g \
    -v -lof -motif -nextProt \
    GRCh37.75 \
    "${SCRIPT_DIR}/data/protocols/ex1.vcf" \
    > "${SCRIPT_DIR}/data/ex1.eff.vcf"

# Run case-control analysis
SnpSift \
    -Xmx${AVAILABLE_MEM_GB}g \
    caseControl -v \
    -tfam "${SCRIPT_DIR}/data/protocols/pedigree.tfam" \
    "${SCRIPT_DIR}/data/ex1.eff.vcf" \
    > "${SCRIPT_DIR}/data/ex1.eff.cc.vcf"

# Filter variants based on case-control and effect impact
cat "${SCRIPT_DIR}/data/ex1.eff.cc.vcf" \
    | SnpSift filter \
    "(Cases[0] = 3) & (Controls[0] = 0) & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" \
    > "${SCRIPT_DIR}/data/ex1.filtered.vcf"

# Download and prepare ClinVar data
wget -O "${SCRIPT_DIR}/data/clinvar_2025.vcf.gz" \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20250521.vcf.gz
gunzip -k "${SCRIPT_DIR}/data/clinvar_2025.vcf.gz"

# Annotate with ClinVar
SnpSift \
    -Xmx${AVAILABLE_MEM_GB}g \
    annotate -v \
    "${SCRIPT_DIR}/data/clinvar_2025.vcf" \
    "${SCRIPT_DIR}/data/ex1.eff.cc.vcf" \
    > "${SCRIPT_DIR}/data/ex1.eff.cc.clinvar.vcf"

# Filter for CFTR variants
cat "${SCRIPT_DIR}/data/ex1.eff.cc.clinvar.vcf" \
    | SnpSift filter \
    "(GENEINFO[*] =~ 'CFTR') & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE')) & (Cases[0] = 3) & (Controls[0] = 0)" \
    > "${SCRIPT_DIR}/output_cf_variant.txt"
