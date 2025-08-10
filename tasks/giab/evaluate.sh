#!/usr/bin/env bash
set -euo pipefail

# Evaluate SNP/indel performance on GIAB HG001 chr20 using hap.py and truvari
# Inputs expected in ./benchmark_data:
#  - HG001_chr20_truth.vcf.gz (and .tbi)
#  - HG001_chr20_highconf.bed
# Reference expected in ./benchmark_ref:
#  - GRCh38_noalt.fa
# Query VCF is auto-detected from ./results unless provided via -q/--query

THREADS=${THREADS:-$(nproc)}
CWD="$(pwd)"

TRUTH_VCF="$CWD/benchmark_data/HG001_chr20_truth.vcf.gz"
REGIONS_BED="$CWD/benchmark_data/HG001_chr20_highconf.bed"
REFERENCE_FA="$CWD/benchmark_ref/GRCh38_noalt.fa"
REFERENCE_SDF_DIR="$CWD/benchmark_ref/GRCh38_noalt.sdf"

OUTDIR="$CWD/evaluation"
HAPPY_DIR="$OUTDIR/happy"
TRUVARI_DIR="$OUTDIR/truvari"

QUERY_VCF=""

usage() {
  echo "Usage: $0 [-q QUERY_VCF]"
  echo "  -q, --query   Path to query VCF.gz to evaluate (if omitted, auto-detect in ./results)"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -q|--query)
      QUERY_VCF="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

echo "Checking required inputs..."
if [[ ! -f "$TRUTH_VCF" ]]; then
  echo "Missing truth VCF: $TRUTH_VCF" >&2
  exit 1
fi
if [[ ! -f "$TRUTH_VCF.tbi" ]]; then
  echo "Missing truth VCF index: $TRUTH_VCF.tbi" >&2
  exit 1
fi
if [[ ! -f "$REGIONS_BED" ]]; then
  echo "Missing confident regions BED: $REGIONS_BED" >&2
  exit 1
fi
if [[ ! -f "$REFERENCE_FA" ]]; then
  echo "Missing reference FASTA: $REFERENCE_FA" >&2
  exit 1
fi

# Auto-detect query VCF if not provided
if [[ -z "$QUERY_VCF" ]]; then
  echo "Auto-detecting query VCF in ./results ..."
  mapfile -t candidates < <(find "$CWD/results" -type f -name "*.vcf.gz" 2>/dev/null | grep -v '\\.g\\.vcf\\.gz$' || true)
  if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "No candidate query VCFs found under ./results. Provide one with -q." >&2
    exit 1
  elif [[ ${#candidates[@]} -eq 1 ]]; then
    QUERY_VCF="${candidates[0]}"
  else
    # Prefer files containing 'pass' or 'filtered', else pick the largest
    mapfile -t preferred < <(printf '%s\n' "${candidates[@]}" | grep -E '/(pass|filtered).*vcf.gz$' || true)
    if [[ ${#preferred[@]} -gt 0 ]]; then
      candidates=("${preferred[@]}")
    fi
    # Pick largest by size
    largest=""; largest_size=0
    for f in "${candidates[@]}"; do
      sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
      if (( sz > largest_size )); then largest_size=$sz; largest="$f"; fi
    done
    QUERY_VCF="$largest"
  fi
  echo "Selected query VCF: $QUERY_VCF"
fi

if [[ ! -f "$QUERY_VCF" ]]; then
  echo "Query VCF not found: $QUERY_VCF" >&2
  exit 1
fi

mkdir -p "$HAPPY_DIR" "$TRUVARI_DIR"

echo "Building RTG SDF (if needed) ..."
if [[ ! -d "$REFERENCE_SDF_DIR" ]]; then
  mamba run -n giab rtg format -o "$REFERENCE_SDF_DIR" "$REFERENCE_FA"
fi

echo "Running hap.py (vcfeval engine) ..."
mamba run -n giab hap.py \
  "$TRUTH_VCF" \
  "$QUERY_VCF" \
  -f "$REGIONS_BED" \
  -r "$REFERENCE_FA" \
  --engine vcfeval \
  --threads "$THREADS" \
  -o "$HAPPY_DIR/happy"

echo "hap.py done. Summary: $HAPPY_DIR/happy.summary.csv"

echo "Running truvari ..."
mamba run -n giab truvari bench \
  -b "$TRUTH_VCF" \
  -c "$QUERY_VCF" \
  -f "$REFERENCE_FA" \
  -o "$TRUVARI_DIR" \
  --includebed "$REGIONS_BED" \
  --passonly \
  --threads "$THREADS"

echo "truvari done. Summary: $TRUVARI_DIR/summary.json"

echo "Finished. Outputs in: $OUTDIR"


