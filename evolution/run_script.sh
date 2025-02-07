#!/bin/bash

# Create individual environments for each tool
conda create -n fastp_env -c bioconda fastp -y
conda create -n fastqc_env -c bioconda fastqc -y
conda create -n multiqc_env -c bioconda multiqc -y
conda create -n spades_env -c bioconda spades -y
conda create -n quast_env -c bioconda quast -y

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
 -I data/evol2_R2.fastq.gz -o trimmed/evol2_R1.fastq.gz -O trimmed/evo2_R2.fastq.gz

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

 quast -o assembl/quast assembly/spades-150/scaffolds.fasta assembly/spades-original/scaffolds.fasta

