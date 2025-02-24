mamba create -n rnaseq
mamba activate rnaseq
mamba install -c bioconda fastqc multiqc

mkdir data
prefetch SRR1278968 SRR1278969 SRR1278970 SRR1278971 SRR1278972 SRR1278973 -O data/ 

mamba install -c bioconda fastqc multiqc

mkdir -p data/processing/0_fasterqdump
for sra_file in data/SRR*/**.sra; do
    fasterq-dump "$sra_file" -O data/processing/0_fasterqdump
done 

mkdir 1_fastqc
fastqc data/processing/0_fasterqdump/*fastq -o data/processing/1_fastqc -t 32
multiqc data/processing/1_fastqc -o data/processing/2_multiqc

mamba install bioconda::trimmomatic
mkdir data/processing/3_trimming

for file in data/processing/0_fasterqdump/*_1.fastq; do
    base=$(basename "$file" _1.fastq) 
    trimmomatic PE -threads 32 \
        data/processing/0_fasterqdump/${base}_1.fastq \
        data/processing/0_fasterqdump/${base}_2.fastq \
        data/processing/2_trimming/${base}_1P.fastq \
        data/processing/2_trimming/${base}_1U.fastq \
        data/processing/2_trimming/${base}_2P.fastq \
        data/processing/2_trimming/${base}_2U.fastq \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done

mkdir data/processing/4_fastqc_trimmed
mkdir data/processing/5_multiqc_trimmed

fastqc data/processing/2_trimming/*fastq -o data/processing/4_fastqc_trimmed -t 32
multiqc data/processing/4_fastqc_trimmed -o data/processing/5_multiqc_trimmed

mkdir data/reference
wget -O data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz http://www.candidagenome.org/download/sequence/C_parapsilosis_CDC317/current/C_parapsilosis_CDC317_current_chromosomes.fasta.gz
wget -O data/reference/C_parapsilosis_CDC317_current_features.gff http://www.candidagenome.org/download/gff/C_parapsilosis_CDC317/C_parapsilosis_CDC317_current_features.gff
gunzip -k data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta.gz

mamba install -c bioconda star gffread
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir data/processing/6_indexing/ \
     --genomeFastaFiles data/reference/C_parapsilosis_CDC317_current_chromosomes.fasta \
     --genomeSAindexNbases 10
gffread data/reference/C_parapsilosis_CDC317_current_features.gff -T \
        -o data/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf


mkdir data/processing/7_mapping
for file in data/processing/0_fasterqdump/*_1.fastq; do
    base=$(basename "$file" _1.fastq)
    STAR --runThreadN 32 --genomeDir data/processing/6_indexing \
        --sjdbGTFfile data/processing/6_indexing/C_parapsilosis_CDC317_current_features.gtf \
        --readFilesIn data/processing/2_trimming/${base}_1P.fastq \
                      data/processing/2_trimming/${base}_2P.fastq \
        --outFileNamePrefix data/processing/7_mapping/${base}_ --outSAMtype BAM \
        SortedByCoordinate --limitBAMsortRAM 100000000000 \
        --quantMode GeneCounts
done

multiqc data/processing/7_mapping -o data/processing/7_mapping

# R scripting starts from here
mamba install conda-forge::r-base
mamba install -c bioconda bioconductor-deseq2 bioconductor-biocparallel bioconductor-clusterprofiler
mamba install -c conda-forge r-ggpubr

mkdir data/processing/8_deseq

# download go annotations 
mkdir data/processing/9_gorich
wget -O data/processing/9_gorich/gene_association.cgd.gz http://www.candidagenome.org/download/go/gene_association.cgd.gz
gunzip -k data/processing/9_gorich/gene_association.cgd.gz

Rscript run_deseq.R