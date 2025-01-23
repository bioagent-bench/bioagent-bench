SRA_FILE="sra_ids.txt"

while read -r SRR_ID; do
  [[ -z "$SRR_ID" ]] && continue
  echo "Prefetching $SRR_ID ..."
  mkdir -p "../${SRR_ID}"
  prefetch "$SRR_ID" -o "../data/${SRR_ID}/"
done < "$SRA_FILE"