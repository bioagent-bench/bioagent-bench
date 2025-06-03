#!/bin/bash
set -e

# Function for logging with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# # Get the absolute path of the script directory
SCRIPT_DIR="/app"
OUTPUT_DIR="${SCRIPT_DIR}/output"

AVAILABLE_MEM_GB=4 

# log "Starting fibrosis analysis pipeline..."

# # Create necessary directories
# log "Creating directories..."
# mkdir -p "${SCRIPT_DIR}/data"
# mkdir -p "${OUTPUT_DIR}"

# # Download reference data
# log "Downloading SNPEff reference data (GRCh37.75)..."
# snpEff download -v GRCh37.75

# # Download and extract source data
# log "Downloading and extracting protocol data..."
# wget -O "${SCRIPT_DIR}/data/protocols.zip" \
#     "http://sourceforge.net/projects/snpeff/files/protocols.zip"
# unzip "${SCRIPT_DIR}/data/protocols.zip" -d "${SCRIPT_DIR}/data/"
# rm "${SCRIPT_DIR}/data/protocols.zip" 

# # Run SNP effect prediction
# log "Running SNP effect prediction..."
# snpEff \
#     -Xmx${AVAILABLE_MEM_GB}g \
#     -v -lof -motif -nextProt \
#     GRCh37.75 \
#     "${SCRIPT_DIR}/data/protocols/ex1.vcf" \
#     > "${OUTPUT_DIR}/ex1.eff.vcf"
# log "SNP effect prediction completed"

# # Run case-control analysis
# log "Running case-control analysis..."
# SnpSift \
#     -Xmx${AVAILABLE_MEM_GB}g \
#     caseControl -v \
#     -tfam "${SCRIPT_DIR}/data/protocols/pedigree.tfam" \
#     "${OUTPUT_DIR}/ex1.eff.vcf" \
#     > "${OUTPUT_DIR}/ex1.eff.cc.vcf"
# log "Case-control analysis completed"

# # Filter variants based on case-control and effect impact
# log "Filtering variants..."
# cat "${OUTPUT_DIR}/ex1.eff.cc.vcf" \
#     | SnpSift filter \
#     "(Cases[0] = 3) & (Controls[0] = 0) & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" \
#     > "${OUTPUT_DIR}/ex1.filtered.vcf"

# Download and prepare ClinVar data
# log "Downloading and preparing ClinVar data..."
# wget -O "${SCRIPT_DIR}/data/clinvar.vcf.gz" \
#     https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
# gunzip -k "${SCRIPT_DIR}/data/clinvar.vcf.gz"

# # Annotate with ClinVar
# log "Annotating with ClinVar..."
# SnpSift \
#     -Xmx${AVAILABLE_MEM_GB}g \
#     annotate -v \
#     "${SCRIPT_DIR}/data/clinvar.vcf" \
#     "${OUTPUT_DIR}/ex1.eff.cc.vcf" \
#     > "${OUTPUT_DIR}/ex1.eff.cc.clinvar.vcf"
# log "ClinVar annotation completed"

# # Filter for CFTR variants
# log "Filtering for CFTR variants..."
# cat "${OUTPUT_DIR}/ex1.eff.cc.clinvar.vcf" \
#     | SnpSift filter \
#     "(GENEINFO[*] =~ 'CFTR') & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE')) & (Cases[0] = 3) & (Controls[0] = 0)" \
#     > "${OUTPUT_DIR}/output_cf_variant.txt"

# log "Pipeline completed successfully!"

# Convert VCF to CSV format
log "Converting VCF to CSV format..."
mkdir -p results
python3 vcf_to_csv.py
log "Conversion completed. Results available in results/cf_variants.csv"