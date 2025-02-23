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
fastqc *fastq.gz -o ./fastqc_initial -t 2
multiqc ./1_fastqc

conda install bioconda::trimmomatic