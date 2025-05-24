#! /bin/bash

mamba create -n single-cell -c bioconda kb-python kallisto
mamba create -n kallisto -c bioconda kallisto
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


# Build kallisto index
mamba activate kallisto
mkdir processed/2_indexing
kallisto index -i processed/2_indexing/Homo_sapiens.GRCh38.cdna.all.index \ 
   data/reference/Homo_sapiens.GRCh38.cdna.all.fa

mamba activate single-cell
kb ref -i processed/2_indexing/Homo_sapiens.GRCh38.cdna.all.index \
   -d human \
   -g processed/2_indexing/t2g.txt

mkdir processed/3_quantification
kb count processed/1_merged/pbmc_1k_v3_S1_R1.fastq.gz \
   processed/1_merged/pbmc_1k_v3_S1_R2.fastq.gz \
   -i processed/2_indexing/Homo_sapiens.GRCh38.cdna.all.index \
   -x 10XV3 \
   -g processed/2_indexing/t2g.txt \
   -t 60 \
   --cellranger \
   -o processed/3_quantification

mamba create -n single-cell-r -c conda-forge r-base r-essentials r-biocmanager==3.20
mamba activate single-cell-r
mamba install -c conda-forge r-tidyverse r-seurat r-matrix r-scales r-rjson r-r2html r-dt

chmod +x processing.R
mkdir -p processed/4_analysis
Rscript processing.R