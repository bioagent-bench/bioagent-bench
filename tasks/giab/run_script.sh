#!/usr/bin/env bash
set -euo pipefail

# Minimal nf-core/sarek germline run + hap.py benchmarking (NA12878 WES)
# Requires: Nextflow + Docker. Downloads test FASTQs, BED, GIAB truth set, and GRCh38 FASTA.

CWD="$(pwd)"
OUTDIR="$CWD/results"
WORKDIR="$CWD/output"
DATADIR="$CWD/data"
BENCHDIR="$CWD/benchmark_ref"
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

mkdir -p "$OUTDIR" "$WORKDIR" "$DATADIR" "$BENCHDIR" "$GIAB_DIR" "$REF_DIR"

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
