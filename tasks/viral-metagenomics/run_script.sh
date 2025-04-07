mamba create -n viral-metagenomics -c bioconda fastp bowtie2 kaiju krona megahit prokka
mamba activate viral-metagenomics

mkdir -p data
wget -O data/Dol1_S19_L001_R1_001.fastq.gz https://osf.io/4x6qs/download
wget -O data/Dol1_S19_L001_R2_001.fastq.gz https://osf.io/z2xed/download

mkdir -p data/reference
wget -O data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/tursiops_truncatus/dna/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz

mkdir -p processing/1_trimming
fastp -i data/Dol1_S19_L001_R1_001.fastq.gz -o processing/1_trimming/Dol1_trimmed_R1.fastq \
 -I data/Dol1_S19_L001_R2_001.fastq.gz -O processing/1_trimming/Dol1_trimmed_R2.fastq \
 --detect_adapter_for_pe --length_required 30 \
 --cut_front --cut_tail --cut_mean_quality 10

gunzip -k data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz
mkdir processing/2_indexing
cd processing/2_indexing
bowtie2-build ../../data/reference/Tursiops_truncatus.turTru1.dna.toplevel.fa Tursiops_truncatus --threads 30
cd ..
mkdir 3_mapping
cd 3_mapping
bowtie2 -x Tursiops_truncatus -1 ../1_trimming/Dol1_trimmed_R1.fastq -2 ../1_trimming/Dol1_trimmed_R2.fastq -S ../3_mapping/dol_map.sam --un-conc ../3_mapping/Dol_reads_unmapped.fastq --threads 30
cd ..

mkdir 4_assembly
cd 4_assembly
megahit -1 ../3_mapping/Dol_reads_unmapped.1.fastq -2 ../3_mapping/Dol_reads_unmapped.2.fastq


cd ..
mkdir 5_classification
cd 5_classification
kaiju-makedb -s viruses

kaiju -t 5_classification/nodes.dmp -f 5_classification/viruses/kaiju_db_viruses.fmi \
 -i 4_assembly/megahit_out/final.contigs.fa -o 5_classification/Dol1_contigs_kaiju.out

kaiju2krona -t 5_classification/nodes.dmp -n 5_classification/names.dmp \
 -i 5_classification/Dol1_contigs_kaiju.out -o 5_classification/Dol1_contigs_kaiju.krona -u
ktImportText -o 5_classification/Dol1_contigs_kaiju.krona.html 5_classification/Dol1_contigs_kaiju.krona

grep '^>' 5_classification/Dol1_contigs_kaiju.out | sed 's/.*len=\([0-9]*\)/\1/' | sort -n | tail -n 1

