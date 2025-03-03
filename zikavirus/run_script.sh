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

mkdir -p processing/2_alignment/feature_count_out
mamba install -c bioconda subread

# unpaired reads
for fq in processing/0_fasterqdump/*.fastq; do
    base=$(basename "$fq" .fastq)
    if [[ ! "$base" =~ _[12]$ ]]; then
        echo "Processing unpaired read: $fq"
        echo "Base name: $base"
        
        echo "Aligning reads from $base to the reference genome"
        mkdir -p processing/2_alignment/unpaired/$base
        
        STAR \
            --genomeDir data/reference/star_index \
            --sjdbGTFfile data/reference/hg38/hg38.knownGene.gtf \
            --runThreadN 32 \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix processing/2_alignment/unpaired/$base/ \
            --readFilesIn $fq \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --outSAMmode Full
            
        suffix="Aligned.out.bam"
        outname="$base.count.txt"
        bam="processing/2_alignment/$base/$suffix"
        
        featureCounts \
            -T 32 \
            -t exon \
            -g gene_id \
            -a data/reference/hg38/hg38.knownGene.gtf \
            -o processing/2_alignment/feature_count_out/$outname \
            $bam
    fi
done

mkdir -p processing/2_alignment/paired

# paired reads
for fq1 in processing/0_fasterqdump/*_1.fastq; do
    basename=$(basename "$fq1" _1.fastq)
    fq2="${fq1/_1.fastq/_2.fastq}"
    
    if [ -f "$fq2" ]; then
        echo "Processing paired-end reads for: $basename"
        
        echo "Aligning reads from $basename to the reference genome"
        mkdir -p processing/2_alignment/paired/$basename
        
        STAR \
            --genomeDir data/reference/star_index \
            --sjdbGTFfile data/reference/hg38/hg38.knownGene.gtf \
            --runThreadN 32 \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix processing/2_alignment/paired/$basename/ \
            --readFilesIn $fq1 $fq2 \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --outSAMmode Full
            
        suffix="Aligned.out.bam"
        outname="$basename.count.txt"
        bam="processing/2_alignment/paired/$basename/$suffix"
        
        featureCounts \
            -T 32 \
            -t exon \
            -g gene_id \
            -a data/reference/hg38/hg38.knownGene.gtf \
            -o processing/2_alignment/feature_count_out/$outname \
            $bam
    else
        echo "Warning: Could not find matching pair file for $fq1"
    fi
done
