#!/bin/bash
set -e

mkdir data
wget -O data/protocols.zip http://sourceforge.net/projects/snpeff/files/protocols.zip
unzip data/protocols.zip -d data/

mamba create -n fibrosis_env -c bioconda snpeff snpsift
mamba activate fibrosis_env
snpEff download -v GRCh37.75
snpEff -Xmx12g -v -lof -motif -nextProt GRCh37.75 data/protocols/ex1.vcf > data/ex1.eff.vcf

awk -F',' 'NR>1 {print "CEPH_1463", $1, ($3=="0" ? "0" : $3), ($4=="0" ? "0" : $4), $5, $6}' family_tree_metadata.csv > data/pedigree.tfam

SnpSift -Xmx100g caseControl -v -tfam data/pedigree.tfam data/ex1.eff.vcf \
  > data/ex1.eff.cc.vcf
