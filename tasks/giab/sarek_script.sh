#!/usr/bin/env bash
set -euo pipefail

# Minimal nf-core/sarek germline variant calling runner
# Assumptions:
# - Paired FASTQs exist at /fastq/HG001_chr20_R1.fq.gz and /fastq/HG001_chr20_R2.fq.gz
# - Reference FASTA exists at ./benchmark_ref/GRCh38_noalt.fa
# - Nextflow and Docker are installed and available on PATH

CWD="$(pwd)"
OUTDIR="$CWD/results"
WORKDIR="$CWD/output"
SAREK_RELEASE="3.5.1"
PROFILE="docker"
SAMPLE_ID="HG001_chr20"
R1="$CWD/data/HG001_chr20_R1.fq.gz"
R2="$CWD/data/HG001_chr20_R2.fq.gz"
THREADS=64

mkdir -p "$OUTDIR" "$WORKDIR"

SAMPLESHEET="$OUTDIR/samplesheet.csv"
echo "patient,status,sample,lane,fastq_1,fastq_2" > "$SAMPLESHEET"
echo "${SAMPLE_ID},0,${SAMPLE_ID},lane_1,${R1},${R2}" >> "$SAMPLESHEET"

set -x
nextflow run nf-core/sarek -r "$SAREK_RELEASE" \
  -profile "$PROFILE" \
  --genome GATK.GRCh38 \
  --input "$SAMPLESHEET" \
  --analysis_type "germline" \
  --tools deepvariant,haplotypecaller \
  --outdir "$OUTDIR" \
  --max_cpus "$THREADS" \
  --max_memory "120.GB" \
  -work-dir "$WORKDIR"
set +x

echo "Done. Sarek results: $OUTDIR"
