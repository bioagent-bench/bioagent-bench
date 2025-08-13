#!/usr/bin/env bash
set -euo pipefail

# Variant calling evaluation script using hap.py benchmarking
# This script evaluates DeepVariant and HaplotypeCaller outputs against GIAB truth set

CWD="$(pwd)"
OUTDIR="$CWD/results"
EVALDIR="$CWD/evaluation"

SAMPLE_ID="A006850052_NA12878_75M"

# Input files (simplified paths)
BED_FILE="$CWD/data/Agilent_v7.chr.bed"
TRUTH_VCF="$CWD/benchmark_ref/GIAB_HG001_GRCh38_v4.2.1/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="$CWD/benchmark_ref/GIAB_HG001_GRCh38_v4.2.1/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
REF_FASTA="$CWD/benchmark_ref/ref_GRCh38/Homo_sapiens_assembly38.fasta"

# Variant calling outputs from Sarek
DEEPVARIANT_VCF="$OUTDIR/variant_calling/deepvariant/$SAMPLE_ID/${SAMPLE_ID}.deepvariant.vcf.gz"
HAPLOTYPECALLER_VCF="$OUTDIR/variant_calling/haplotypecaller/$SAMPLE_ID/${SAMPLE_ID}.haplotypecaller.filtered.vcf.gz"

# Environment name for hap.py
HAP_ENV="hap_env"

mkdir -p "$EVALDIR"

echo "== Checking required input files =="
for file in "$TRUTH_VCF" "$TRUTH_BED" "$REF_FASTA" "$BED_FILE" "$DEEPVARIANT_VCF" "$HAPLOTYPECALLER_VCF"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        echo "Please run run_script.sh first to generate variant calling outputs"
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

echo "== Evaluation Complete =="
echo "Summary report: $SUMMARY_REPORT"
echo ""
echo "Key output files:"
echo "- Summary report: $SUMMARY_REPORT"
echo "- DeepVariant evaluation: $EVALDIR/${SAMPLE_ID}_deepvariant.*"
echo "- HaplotypeCaller evaluation: $EVALDIR/${SAMPLE_ID}_haplotypecaller.*"
echo ""
echo "To view the summary:"
echo "cat $SUMMARY_REPORT"
echo ""
echo "To examine detailed metrics, check the .summary.csv files in $EVALDIR/"
