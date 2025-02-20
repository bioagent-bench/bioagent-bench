#!/bin/bash
set -e

mkdir data
# This is the source data from the tutorial
# TODO: Replace this download link with a real upload folder
wget -O data/fibrosis.tar.gz DL_LINK 
tar -xzf data/fibrosis.tar.gz -C data/

mamba create -n fibrosis_env -c bioconda snpeff snpsift
mamba activate fibrosis_env
snpEff download -v GRCh37.75
snpEff -Xmx12g -v -lof -motif -nextProt GRCh37.75 data/protocols/ex1.vcf > data/ex1.eff.vcf

# this converts the csv metadata to tfam format
awk -F',' 'NR>1 {print "CEPH_1463", $1, ($3=="0" ? "0" : $3), ($4=="0" ? "0" : $4), $5, $6}' family_tree_metadata.csv > data/pedigree.tfam

SnpSift -Xmx100g caseControl -v -tfam data/pedigree.tfam data/ex1.eff.vcf \
  > data/ex1.eff.cc.vcf

cat data/ex1.eff.cc.vcf \
| SnpSift filter \
"(Cases[0] = 3) & (Controls[0] = 0) & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" \
> data/ex1.filtered.vcf

wget -O data/clinvar_2025.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20250217.vcf.gz
gunzip -k data/clinvar_2025.vcf.gz

SnpSift -Xmx100g annotate -v data/clinvar_2025.vcf data/ex1.eff.cc.vcf \
> data/ex1.eff.cc.clinvar.vcf

cat data/ex1.eff.cc.clinvar.vcf \
| SnpSift filter \
    "(GENEINFO[*] =~ 'CFTR') & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE')) & (Cases[0] = 3) & (Controls[0] = 0)" \
> output_cf_variant.txt

