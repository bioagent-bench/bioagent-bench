mkdir data
wget -O data/metagenomics.zip https://zenodo.org/records/7010950/files/dc_workshop.zip
wget -O data/MGRAST_MetaData_JP.xlsx https://zenodo.org/records/7010950/files/MGRAST_MetaData_JP.xlsx

unzip data/metagenomics.zip -d data/

mamba create -n metagenomics
mamba activate metagenomics

mamba install -c bioconda fastqc multiqc

mkdir -p data/processing/0_fastqc
fastqc data/seq/*fastq.gz -o data/processing/0_fastqc -t 32

mamba install -c bioconda trimmomatic
mkdir data/processing/1_trimmed

for file in data/seq/*.fastq.gz; do
    gunzip -k $file
done

for file in data/seq/*_R1.fastq; do
    base=$(basename "$file" _R1.fastq) 
    trimmomatic PE -threads 32 \
        data/seq/${base}_R1.fastq \
        data/seq/${base}_R2.fastq \
        data/processing/1_trimmed/${base}_1P.fastq \
        data/processing/1_trimmed/${base}_1U.fastq \
        data/processing/1_trimmed/${base}_2P.fastq \
        data/processing/1_trimmed/${base}_2U.fastq \
        SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:data/seq/TruSeq3-PE.fa:2:40:15
done

mkdir -p data/processing/2_fastqc
fastqc data/processing/1_trimmed/*.fastq -o data/processing/2_fastqc -t 32
multiqc data/processing/2_fastqc -o data/processing/2_fastqc

mkdir data/processing/3_assembly
mamba install -c bioconda spades

metaspades.py -1 data/processing/1_trimmed/JC1A_1P.fastq \
              -2 data/processing/1_trimmed/JC1A_2P.fastq \
              -o data/processing/3_assembly/assembly_JC1A

metaspades.py -1 data/processing/1_trimmed/JP4D_1P.fastq \
              -2 data/processing/1_trimmed/JP4D_2P.fastq \
              -o data/processing/3_assembly/assembly_JP4D

# TODO: add metaquast for spades QC

mamba install -c bioconda maxbin2
mkdir data/processing/4_binning
# this will not work as the sample size is too small but thats expected
run_MaxBin.pl -thread 32 \
              -contig data/processing/3_assembly/assembly_JC1A/contigs.fasta \
              -reads data/processing/1_trimmed/JC1A_1P.fastq \
              -reads2 data/processing/1_trimmed/JC1A_2P.fastq \
              -out data/processing/4_binning/

run_MaxBin.pl -thread 32 \
              -contig data/processing/3_assembly/assembly_JP4D/contigs.fasta \
              -reads data/processing/1_trimmed/JP4D_1P.fastq \
              -reads2 data/processing/1_trimmed/JP4D_2P.fastq \
              -out data/processing/4_binning/

    