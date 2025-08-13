#!/usr/bin/env bash
set -euo pipefail

# Minimal nf-core/sarek germline run + hap.py benchmarking (NA12878 WES)
# Requires: Nextflow + Docker. Downloads test FASTQs, BED, GIAB truth set, and GRCh38 FASTA.
# Includes automated evaluation using hap.py after variant calling completes.

CWD="$(pwd)"
OUTDIR="$CWD/results"
WORKDIR="$CWD/output"
DATADIR="$CWD/data"
BENCHDIR="$CWD/benchmark_ref"
EVALDIR="$CWD/evaluation"
SAREK_RELEASE="3.5.1"
PROFILE="docker"

SAMPLE_ID="A006850052_NA12878_75M"
R1="$DATADIR/A006850052_NA12878_75M_R1.fq.gz"
R2="$DATADIR/A006850052_NA12878_75M_R2.fq.gz"
BED_FILE="$DATADIR/Agilent_v7.bed"         # panel/exome targets (restricts Sarek + benchmarking)
THREADS=${THREADS:-64}

# -----------------------------
# Data URLs (Zenodo NA12878 Agilent 75M WES)
# -----------------------------
R1_URL="https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R1.fq.gz?download=1"
R2_URL="https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R2.fq.gz?download=1"
BED_URL="https://zenodo.org/records/6513789/files/Agilent_v7.bed?download=1"

# -----------------------------
# GIAB truth set (HG001 / NA12878, GRCh38 v4.2.1)
# -----------------------------
GIAB_DIR="$BENCHDIR/GIAB_HG001_GRCh38_v4.2.1"
TRUTH_VCF="$GIAB_DIR/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="$GIAB_DIR/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
TRUTH_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38"
# -----------------------------
# Reference FASTA (GRCh38)
# If you already have a matching GRCh38 FASTA, set REF_FASTA to it and skip download.
# -----------------------------
REF_DIR="$BENCHDIR/ref_GRCh38"
REF_FASTA="${REF_FASTA:-$REF_DIR/Homo_sapiens_assembly38.fasta}"
REF_FASTA_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
REF_DICT_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"  # dict is small
# (fai will be built locally; we also work fine if you already have .fai / bwa indexes)

mkdir -p "$OUTDIR" "$WORKDIR" "$DATADIR" "$BENCHDIR" "$GIAB_DIR" "$REF_DIR" "$EVALDIR"

echo "== Checking input WES dataset =="
if [[ ! -f "$R1" ]]; then wget -O "$R1" "$R1_URL"; fi
if [[ ! -f "$R2" ]]; then wget -O "$R2" "$R2_URL"; fi
if [[ ! -f "$BED_FILE" ]]; then wget -O "$BED_FILE" "$BED_URL"; fi

BED_FILE_CHR="$DATADIR/Agilent_v7.chr.bed"
awk 'BEGIN{OFS="\t"}
     /^track|^browser|^#/ {print; next}
     { if ($1 !~ /^chr/) $1="chr"$1; print }' "$BED_FILE" > "$BED_FILE_CHR"
BED_FILE="$BED_FILE_CHR"

echo "== Fetching GIAB truth set (HG001 GRCh38 v4.2.1) =="
if [[ ! -f "$TRUTH_VCF" ]]; then wget -O "$TRUTH_VCF" "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"; fi
if [[ ! -f "${TRUTH_VCF}.tbi" ]]; then wget -O "${TRUTH_VCF}.tbi" "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"; fi
if [[ ! -f "$TRUTH_BED" ]]; then wget -O "$TRUTH_BED" "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"; fi

echo "== Preparing GRCh38 reference for hap.py =="
if [[ ! -f "$REF_FASTA" ]]; then
  # Download the GRCh38 FASTA (GATK bundle mirror on GCS) - it's already uncompressed
  wget -O "$REF_FASTA" "$REF_FASTA_URL"
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then samtools faidx "$REF_FASTA"; fi
# Optional: dictionary (not strictly required by hap.py, but handy)
if [[ ! -f "${REF_DIR}/Homo_sapiens_assembly38.dict" ]]; then
  wget -O "${REF_DIR}/Homo_sapiens_assembly38.dict" "$REF_DICT_URL" || true
fi

echo "== Build Sarek samplesheet =="
SAMPLESHEET="$OUTDIR/samplesheet.csv"
echo "patient,status,sample,lane,fastq_1,fastq_2" > "$SAMPLESHEET"
echo "${SAMPLE_ID},0,${SAMPLE_ID},lane_1,${R1},${R2}" >> "$SAMPLESHEET"

echo "== Run Sarek (germline, WES, intervals) =="
set -x
nextflow run nf-core/sarek -r "$SAREK_RELEASE" \
  -profile "$PROFILE" \
  --fasta "$REF_FASTA" \
  --input "$SAMPLESHEET" \
  --analysis_type germline \
  --tools deepvariant,haplotypecaller \
  --wes \
  --intervals "$BED_FILE" \
  --outdir "$OUTDIR" \
  --max_cpus "$THREADS" \
  --max_memory "120.GB" \
  -work-dir "$WORKDIR"
set +x

echo "Sarek results:        $OUTDIR"

echo ""
echo "==============================================="
echo "Starting variant calling evaluation with hap.py"
echo "==============================================="

# Define paths to variant calling outputs from Sarek
DEEPVARIANT_VCF="$OUTDIR/variant_calling/deepvariant/$SAMPLE_ID/${SAMPLE_ID}.deepvariant.vcf.gz"
HAPLOTYPECALLER_VCF="$OUTDIR/variant_calling/haplotypecaller/$SAMPLE_ID/${SAMPLE_ID}.haplotypecaller.filtered.vcf.gz"

# Environment name for hap.py
HAP_ENV="hap_env"

echo "== Checking required input files for evaluation =="
for file in "$TRUTH_VCF" "$TRUTH_BED" "$REF_FASTA" "$BED_FILE" "$DEEPVARIANT_VCF" "$HAPLOTYPECALLER_VCF"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        echo "Variant calling may have failed or files are in unexpected locations"
        exit 1
    fi
done

echo "== Setting up hap.py environment =="
# Check if mamba is available, fallback to conda
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "ERROR: Neither mamba nor conda found. Please install mamba or conda first."
    exit 1
fi

# Check if environment exists, create if not
if ! $CONDA_CMD env list | grep -q "^${HAP_ENV}\s"; then
    echo "Creating hap.py environment..."
    $CONDA_CMD create -n "$HAP_ENV" -c bioconda -c conda-forge hap.py rtg-tools -y
else
    echo "hap.py environment already exists"
fi

# Function to run hap.py evaluation
run_happy_eval() {
    local vcf_file="$1"
    local caller_name="$2"
    local output_prefix="$EVALDIR/${SAMPLE_ID}_${caller_name}"
    
    echo "== Running hap.py evaluation for $caller_name =="
    echo "Input VCF: $vcf_file"
    echo "Output prefix: $output_prefix"
    
    $CONDA_CMD run -n "$HAP_ENV" hap.py \
        "$TRUTH_VCF" \
        "$vcf_file" \
        -f "$TRUTH_BED" \
        -T "$BED_FILE" \
        -r "$REF_FASTA" \
        -o "$output_prefix" \
        --pass-only \
        --verbose
    
    echo "Completed evaluation for $caller_name"
    echo "Results saved to: ${output_prefix}.*"
}

# Run evaluations for both variant callers
run_happy_eval "$DEEPVARIANT_VCF" "deepvariant"
run_happy_eval "$HAPLOTYPECALLER_VCF" "haplotypecaller"

echo "== Generating summary report =="
SUMMARY_REPORT="$EVALDIR/evaluation_summary.txt"
{
    echo "========================================="
    echo "Variant Calling Evaluation Summary"
    echo "Sample: $SAMPLE_ID"
    echo "Date: $(date)"
    echo "========================================="
    echo ""
    
    for caller in deepvariant haplotypecaller; do
        summary_file="$EVALDIR/${SAMPLE_ID}_${caller}.summary.csv"
        if [[ -f "$summary_file" ]]; then
            echo "--- $caller Results ---"
            echo "Summary file: $summary_file"
            # Extract key metrics from summary
            if [[ -f "$summary_file" ]]; then
                echo "Key metrics:"
                head -n 20 "$summary_file" | grep -E "(Type|TRUTH|QUERY|Recall|Precision|F1_Score)" || echo "Summary format may differ"
            fi
            echo ""
        else
            echo "--- $caller Results ---"
            echo "WARNING: Summary file not found: $summary_file"
            echo ""
        fi
    done
    
    echo "Detailed results available in:"
    echo "- $EVALDIR/"
    echo ""
    echo "Files generated:"
    ls -la "$EVALDIR"/ || true
    
} > "$SUMMARY_REPORT"

echo "== Pipeline Complete =="
echo ""
echo "========================================="
echo "All results:"
echo "========================================="
echo "- Sarek results:      $OUTDIR"
echo "- Evaluation results: $EVALDIR"
echo "- Summary report:     $SUMMARY_REPORT"
echo ""
echo "Key evaluation files:"
echo "- DeepVariant evaluation:     $EVALDIR/${SAMPLE_ID}_deepvariant.*"
echo "- HaplotypeCaller evaluation: $EVALDIR/${SAMPLE_ID}_haplotypecaller.*"
echo ""
echo "To view the evaluation summary:"
echo "cat $SUMMARY_REPORT"
echo ""
echo "To examine detailed metrics, check the .summary.csv files in $EVALDIR/"
