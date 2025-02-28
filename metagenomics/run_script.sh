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

mamba activate metagenomics
mamba install -c bioconda kraken2
mkdir data/processing/4_taxonomy
mkdir data/kraken_db

wget -v -O data/k2_standard_16gb_20241228.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz
tar -xvf data/k2_standard_16gb_20241228.tar.gz -C data/kraken_db

kraken2 --db data/kraken_db --threads 32 \
        data/processing/3_assembly/assembly_JC1A/contigs.fasta \
        --output "data/processing/4_taxonomy/JC1A.kraken" \
        --report "data/processing/4_taxonomy/JC1A.report"

kraken2 --db data/kraken_db --threads 32 \
        data/processing/3_assembly/assembly_JP4D/contigs.fasta \
        --output "data/processing/4_taxonomy/JP4D.kraken" \
        --report "data/processing/4_taxonomy/JP4D.report"

# optional  generate a visualization report with krona
wget -v -O data/taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

mamba install -c bioconda krona
mkdir data/processing/5_krona
cut -f2,3 data/processing/4_taxonomy/JP4D.kraken > data/processing/5_krona/JP4D_krona.input
cut -f2,3 data/processing/4_taxonomy/JC1A.kraken > data/processing/5_krona/JC1A_krona.input
ktImportTaxonomy data/processing/5_krona/JP4D_krona.input -o data/processing/5_krona/JP4D_krona.out.html
ktImportTaxonomy data/processing/5_krona/JC1A_krona.input -o data/processing/5_krona/JC1A_krona.out.html