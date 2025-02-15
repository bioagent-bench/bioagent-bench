#!/bin/bash

# Create individual environments for each tool
# Sometimes these environments get spastic and this is easier
conda create -n fastp_env -c bioconda fastp -y
conda create -n fastqc_env -c bioconda fastqc -y
conda create -n multiqc_env -c bioconda multiqc -y
conda create -n spades_env -c bioconda spades -y
conda create -n quast_env
conda create -y -n mapping_env samtools bwa qualimap r-base

# create the folders and download the data
mkdir analysis
cd analysis
wget -O data.tar.gz https://osf.io/2jc4a/download
tar -xvzf data.tar.gz

# let conda solve the whole environment
conda env create -f ../environment.yml

# QUALITY CONTROL
# trim the adapters
conda activate fastp_env
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 60 \
 --html trimmed/anc.fastp.html --json trimmed/anc.fastp.json -i data/anc_R1.fastq.gz \
 -I data/anc_R2.fastq.gz -o trimmed/anc_R1.fastq.gz -O trimmed/anc_R2.fastq.gz

 fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 60 \
 --html trimmed/evol1.fastp.html --json trimmed/evol1.fastp.json -i data/evol1_R1.fastq.gz \
 -I data/evol1_R2.fastq.gz -o trimmed/evol1_R1.fastq.gz -O trimmed/evol1_R2.fastq.gz

 fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 60 \
 --html trimmed/evol2.fastp.html --json trimmed/evol2.fastp.json -i data/evol2_R1.fastq.gz \
 -I data/evol2_R2.fastq.gz -o trimmed/evol2_R1.fastq.gz -O trimmed/evol2_R2.fastq.gz

 # run fastqc
 conda activate fastqc_env
 fastqc -o trimmed-fastqc trimmed/*.fastq.gz

 # (optional - visualization) visualize the results using multiqc
 # !you need to install multiqc
 conda activate multiqc_env
 multiqc trimmed-fastqc trimmed

 # GENOME ASSEMBLY
 conda activate spades_env
 mkdir assembly
 spades.py -o assembly/spades-150/ --careful -1 trimmed/anc_R1.fastq.gz -2 trimmed/anc_R2.fastq.gz
 spades.py -o assembly/spades-original/ --careful -1 data/anc_R1.fastq.gz -2 data/anc_R2.fastq.gz

# Check the assembly quality
# Programing gods forgive this black magic
 conda activate quast_env
 conda install pip
 pip install setuptools
 pip install quast
 quast -o assembly/quast assembly/spades-150/scaffolds.fasta assembly/spades-original/scaffolds.fasta

# READ MAPPING
conda create -y -n mapping_env samtools bwa qualimap r-base
conda activate mapping_env
mkdir mappings
cp assembly/spades-150/scaffolds.fasta mappings/scaffolds.fasta
bwa index mappings/scaffolds.fasta
bwa mem mappings/scaffolds.fasta trimmed/evol1_R1.fastq.gz trimmed/evol1_R2.fastq.gz > mappings/evol1.sam
bwa mem mappings/scaffolds.fasta trimmed/evol2_R1.fastq.gz trimmed/evol2_R2.fastq.gz > mappings/evol2.sam

# POST PROCESS MAPPING
# Evol1 line
samtools sort -n -O sam mappings/evol1.sam | samtools fixmate -m -O bam - mappings/evol1.fixmate.bam
# sort
samtools sort -O bam -o mappings/evol1.sorted.bam mappings/evol1.fixmate.bam
# mark duplicates
samtools markdup -r -S mappings/evol1.sorted.bam mappings/evol1.sorted.dedup.bam
# extract q20 mappers
samtools view -h -b -q 20 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.dedup.q20.bam
# extract unmapped
samtools view -b -f 4 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.unmapped.bam
# covert to fastq
samtools fastq -1 mappings/evol1.sorted.unmapped.R1.fastq.gz -2 mappings/evol1.sorted.unmapped.R2.fastq.gz mappings/evol1.sorted.unmapped.bam
# delete unnecessary files
rm mappings/evol1.sam
rm mappings/evol1.fixmate.bam
rm mappings/evol1.sorted.bam
rm mappings/evol1.sorted.dedup.bam
rm mappings/evol1.sorted.unmapped.bam

# Evol2 line
samtools sort -n -O sam mappings/evol2.sam | samtools fixmate -m -O bam - mappings/evol2.fixmate.bam
samtools sort -O bam -o mappings/evol2.sorted.bam mappings/evol2.fixmate.bam
samtools markdup -r -S mappings/evol2.sorted.bam mappings/evol2.sorted.dedup.bam
samtools view -h -b -q 20 mappings/evol2.sorted.dedup.bam > mappings/evol2.sorted.dedup.q20.bam
samtools view -b -f 4 mappings/evol2.sorted.dedup.bam > mappings/evol2.sorted.unmapped.bam
samtools fastq -1 mappings/evol2.sorted.unmapped.R1.fastq.gz -2 mappings/evol2.sorted.unmapped.R2.fastq.gz mappings/evol2.sorted.unmapped.bam

rm mappings/evol2.sam
rm mappings/evol2.fixmate.bam
rm mappings/evol2.sorted.bam
rm mappings/evol2.sorted.dedup.bam
rm mappings/evol2.sorted.unmapped.bam

# TAXONOMIC CLASSIFICATION
$ conda create -y -n kraken_env kraken2 bracken
$ conda activate kraken_env
mkdir kraken-db
wget -O kraken-db/k2_standard_16gb_20241228.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz
kraken2 --use-names --threads 32 --db kraken-db/ --fastq-input --report evol1 --gzip-compressed --paired mappings/evol1.sorted.unmapped.R1.fastq.gz mappings/evol1.sorted.unmapped.R2.fastq.gz > evol1.kraken
bracken -d kraken-db -i evol1.kraken -l S -o evol1.bracken

# VARIANT CALLING
conda create -y -n variant_env samtools bamtools freebayes bedtools vcflib rtg-tools bcftools matplotlib
mkdir variants
samtools faidx assembly/scaffolds.fasta
bamtools index -in mappings/evol1.sorted.dedup.q20.bam
freebayes -p 1 -f assembly/spades-150/scaffolds.fasta mappings/evol1.sorted.dedup.q20.bam > variants/evol1.freebayes.vcf
bgzip variants/evol1.freebayes.vcf
tabix -p vcf variants/evol1.freebayes.vcf.gz
rtg vcfstats variants/evol1.freebayes.vcf.gz