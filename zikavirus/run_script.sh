#!/bin/bash
set -e

mamba create -n zika
mamba activate zika
mamba install -c bioconda fastqc multiqc
mkdir data

prefetch SRX1605080 SRX1605079 SRX1605078 SRX1605077 SRX1602857 SRX1602856 SRX1602855 SRX1602854 -O data --order kart

mkdir -p processing/0_fasterqdump
for sra_file in data/SRR*/**.sra; do
    fasterq-dump "$sra_file" -O processing/0_fasterqdump
done 

mkdir processing/1_fastqc
fastqc processing/0_fasterqdump/*fastq -o processing/1_fastqc -t 32

multiqc processing/1_fastqc -o processing/1_fastqc

mkdir data/reference/hg38
wget -O data/reference/hg38/hg38.fa.gz \
    'https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz'

wget -O data/reference/hg38/hg38.knownGene.gtf.gz \
     'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz'

mamba install -c bioconda star
mkdir data/reference/star_index

gunzip -c data/reference/hg38/hg38.fa.gz > data/reference/hg38/hg38.fa
gunzip -c data/reference/hg38/hg38.knownGene.gtf.gz > data/reference/hg38/hg38.knownGene.gtf

STAR --runThreadN 32 \
    --runMode genomeGenerate \
    --genomeDir data/reference/star_index \
    --genomeFastaFiles data/reference/hg38/hg38.fa \
    --sjdbGTFfile data/reference/hg38/hg38.knownGene.gtf \
    --sjdbOverhang 100
