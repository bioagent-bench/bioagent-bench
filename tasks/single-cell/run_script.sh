#! /bin/bash

mamba create -n single-cell -c bioconda kb-python
mamba activate single-cell

mkdir -p data
mkdir -p processed/1_merged
mkdir -p data/reference

# Download the fastq files
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar -O data/pbmc_1k_v3_fastqs.tar

tar -xf data/pbmc_1k_v3_fastqs.tar -C data/

# Download reference data
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O data/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip -k data/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Merge lane files
pigz -dc data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
     data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz |
  pigz > processed/1_merged/pbmc_1k_v3_S1_R1.fastq.gz

pigz -dc data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
     data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz |
  pigz > processed/1_merged/pbmc_1k_v3_S1_R2.fastq.gz

pigz -dc data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz \
     data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz |
  pigz > processed/1_merged/pbmc_1k_v3_S1_I1.fastq.gz
