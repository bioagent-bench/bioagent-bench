#!/bin/bash
set -e

# Set dynamic thread count
THREADS=$(nproc)

env_exists() {
    local env_name=$1
    mamba env list | grep -q "^${env_name} "
}

# Sometimes these environments get spastic and this is easier
if ! env_exists "fastp_env"; then
    mamba create -n fastp_env bioconda::fastp -y
fi

if ! env_exists "fastqc_env"; then
    mamba create -n fastqc_env bioconda::fastqc -y
fi

if ! env_exists "multiqc_env"; then
    mamba create -n multiqc_env bioconda::multiqc -y
fi

if ! env_exists "spades_env"; then
    mamba create -n spades_env bioconda::spades -y
fi

if ! env_exists "quast_env"; then
    mamba create -n quast_env -y
fi

if ! env_exists "mapping_env"; then
    mamba create -y -n mapping_env -c bioconda samtools bwa
fi

if ! env_exists "variant_env"; then
    mamba create -y -n variant_env -c bioconda samtools bamtools freebayes
fi

if ! env_exists "compleasm_env"; then
    mamba create -y -n compleasm_env bioconda::compleasm
fi

if ! env_exists "prokka_env"; then
    mamba create -y -n prokka_env bioconda::prokka --no-plugins
fi

if ! env_exists "voi_env"; then
    mamba create -n voi_env -c bioconda snpeff genometools-genometools --no-plugins
fi

# Create directory structure
mkdir -p data
mkdir -p /outputs/trimmed
mkdir -p /outputs/trimmed-fastqc
mkdir -p /outputs/assembly
mkdir -p /outputs/mappings
mkdir -p /outputs/variants
mkdir -p /outputs/annotation
mkdir -p /outputs/voi
mkdir -p /results

# Download the data to data folder
cd data
wget -O data.tar.gz https://osf.io/2jc4a/download
tar -xvzf data.tar.gz
cd ..

# QUALITY CONTROL
# trim the adapters
mamba activate fastp_env
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html /outputs/trimmed/anc.fastp.html --json /outputs/trimmed/anc.fastp.json -i data/data/anc_R1.fastq.gz \
 -I data/data/anc_R2.fastq.gz -o /outputs/trimmed/anc_R1.fastq.gz -O /outputs/trimmed/anc_R2.fastq.gz

fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html /outputs/trimmed/evol1.fastp.html --json /outputs/trimmed/evol1.fastp.json -i data/data/evol1_R1.fastq.gz \
 -I data/data/evol1_R2.fastq.gz -o /outputs/trimmed/evol1_R1.fastq.gz -O /outputs/trimmed/evol1_R2.fastq.gz

fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html /outputs/trimmed/evol2.fastp.html --json /outputs/trimmed/evol2.fastp.json -i data/data/evol2_R1.fastq.gz \
 -I data/data/evol2_R2.fastq.gz -o /outputs/trimmed/evol2_R1.fastq.gz -O /outputs/trimmed/evol2_R2.fastq.gz

 # run fastqc
mamba activate fastqc_env
fastqc -t $THREADS -o /outputs/trimmed-fastqc /outputs/trimmed/*.fastq.gz

# (optional - visualization) visualize the results using multiqc
# !you need to install multiqc
mamba activate multiqc_env
multiqc /outputs/trimmed-fastqc /outputs/trimmed -o /outputs

# GENOME ASSEMBLY
mamba activate spades_env
spades.py -o /outputs/assembly --careful -t $THREADS -1 /outputs/trimmed/anc_R1.fastq.gz -2 /outputs/trimmed/anc_R2.fastq.gz

# Check the assembly quality
# Programing gods forgive this black magic
mamba activate quast_env
mamba install pip
pip install setuptools
pip install quast
quast -t $THREADS -o /outputs/assembly/quast /outputs/assembly/scaffolds.fasta /outputs/assembly/spades-original/scaffolds.fasta

# READ MAPPING
mamba activate mapping_env
cp /outputs/assembly/scaffolds.fasta /outputs/mappings/scaffolds.fasta
bwa index /outputs/mappings/scaffolds.fasta
bwa mem -t $THREADS /outputs/mappings/scaffolds.fasta /outputs/trimmed/evol1_R1.fastq.gz /outputs/trimmed/evol1_R2.fastq.gz > /outputs/mappings/evol1.sam
bwa mem -t $THREADS /outputs/mappings/scaffolds.fasta /outputs/trimmed/evol2_R1.fastq.gz /outputs/trimmed/evol2_R2.fastq.gz > /outputs/mappings/evol2.sam

# POST PROCESS MAPPING
# Evol1 line
samtools sort -@ $THREADS -n -O sam /outputs/mappings/evol1.sam | samtools fixmate -m -O bam - /outputs/mappings/evol1.fixmate.bam
samtools sort -@ $THREADS -O bam -o /outputs/mappings/evol1.sorted.bam /outputs/mappings/evol1.fixmate.bam
samtools markdup -@ $THREADS -r -S /outputs/mappings/evol1.sorted.bam /outputs/mappings/evol1.sorted.dedup.bam
samtools view -@ $THREADS -h -b -q 20 /outputs/mappings/evol1.sorted.dedup.bam > /outputs/mappings/evol1.sorted.dedup.q20.bam
samtools view -@ $THREADS -b -f 4 /outputs/mappings/evol1.sorted.dedup.bam > /outputs/mappings/evol1.sorted.unmapped.bam
samtools fastq -@ $THREADS -1 /outputs/mappings/evol1.sorted.unmapped.R1.fastq.gz -2 /outputs/mappings/evol1.sorted.unmapped.R2.fastq.gz /outputs/mappings/evol1.sorted.unmapped.bam
rm /outputs/mappings/evol1.sam
rm /outputs/mappings/evol1.fixmate.bam
rm /outputs/mappings/evol1.sorted.bam
rm /outputs/mappings/evol1.sorted.dedup.bam
rm /outputs/mappings/evol1.sorted.unmapped.bam

# Evol2 line
samtools sort -@ $THREADS -n -O sam /outputs/mappings/evol2.sam | samtools fixmate -m -O bam - /outputs/mappings/evol2.fixmate.bam
samtools sort -@ $THREADS -O bam -o /outputs/mappings/evol2.sorted.bam /outputs/mappings/evol2.fixmate.bam
samtools markdup -@ $THREADS -r -S /outputs/mappings/evol2.sorted.bam /outputs/mappings/evol2.sorted.dedup.bam
samtools view -@ $THREADS -h -b -q 20 /outputs/mappings/evol2.sorted.dedup.bam > /outputs/mappings/evol2.sorted.dedup.q20.bam
samtools view -@ $THREADS -b -f 4 /outputs/mappings/evol2.sorted.dedup.bam > /outputs/mappings/evol2.sorted.unmapped.bam
samtools fastq -@ $THREADS -1 /outputs/mappings/evol2.sorted.unmapped.R1.fastq.gz -2 /outputs/mappings/evol2.sorted.unmapped.R2.fastq.gz /outputs/mappings/evol2.sorted.unmapped.bam

rm /outputs/mappings/evol2.sam
rm /outputs/mappings/evol2.fixmate.bam
rm /outputs/mappings/evol2.sorted.bam
rm /outputs/mappings/evol2.sorted.dedup.bam
rm /outputs/mappings/evol2.sorted.unmapped.bam

# VARIANT CALLING
mamba activate variant_env
samtools faidx /outputs/assembly/scaffolds.fasta
bamtools index -in /outputs/mappings/evol1.sorted.dedup.q20.bam

# evol 1
freebayes -p 1 -f /outputs/assembly/scaffolds.fasta /outputs/mappings/evol1.sorted.dedup.q20.bam > /outputs/variants/evol1.freebayes.vcf
bgzip /outputs/variants/evol1.freebayes.vcf
tabix -p vcf /outputs/variants/evol1.freebayes.vcf.gz
zcat /outputs/variants/evol1.freebayes.vcf.gz | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" | bgzip > /outputs/variants/evol1.freebayes.filtered.vcf.gz
tabix -p vcf /outputs/variants/evol1.freebayes.filtered.vcf.gz

# evol 2
bamtools index -in /outputs/mappings/evol2.sorted.dedup.q20.bam
freebayes -p 1 -f /outputs/assembly/scaffolds.fasta /outputs/mappings/evol2.sorted.dedup.q20.bam > /outputs/variants/evol2.freebayes.vcf
bgzip /outputs/variants/evol2.freebayes.vcf
tabix -p vcf /outputs/variants/evol2.freebayes.vcf.gz
zcat /outputs/variants/evol2.freebayes.vcf.gz | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" | bgzip > /outputs/variants/evol2.freebayes.filtered.vcf.gz
tabix -p vcf /outputs/variants/evol2.freebayes.filtered.vcf.gz

# GENOME ANNOTATION
mamba activate compleasm_env
compleasm run -a /outputs/assembly/scaffolds.fasta -o /outputs/annotation -l bacteria_odb10 -m busco -t $THREADS

mamba activate prokka_env
prokka --cpus $THREADS --kingdom Bacteria --genus Escherichia --species coli --outdir /outputs/annotation/prokka /outputs/assembly/spades-150/scaffolds.fasta

mamba activate voi_env
mkdir -p /outputs/voi/data/mygenome
cp /outputs/assembly/scaffolds.fasta /outputs/voi/data/mygenome/sequences.fa
gzip /outputs/voi/data/mygenome/sequences.fa
cp /outputs/annotation/prokka/*.gff /outputs/voi/data/mygenome/genes.gff
gzip /outputs/voi/data/mygenome/genes.gff

SNPEFF_CONFIG=$(find ~ -name snpEff.config | head -n 1)
cp "$SNPEFF_CONFIG" /outputs/voi/snpEff.config
sed -i 's|^data.dir.*|data.dir = data/|' /outputs/voi/snpEff.config

awk '/# Databases & Genomes/ {print; in_section=1; next} 
     in_section && /^#[-]+$/ {print; print "mygenome.genome : EColiMut"; in_section=0; next} 
     {print}' /outputs/voi/snpEff.config > /outputs/voi/snpEff.config.tmp && mv /outputs/voi/snpEff.config.tmp /outputs/voi/snpEff.config

awk '!seen["mygenome.codonTable"] && !seen["mygenome.checkProtein"] && !seen["mygenome.checkCds"] {print} 
     /# Databases & Genomes/ {seen["mygenome.codonTable"]=1; seen["mygenome.checkProtein"]=1; seen["mygenome.checkCds"]=1} 
     END {print "\nmygenome.codonTable : Bacterial\nmygenome.checkProtein : false\nmygenome.checkCds : false"}' \
     /outputs/voi/snpEff.config > /outputs/voi/snpEff.config.tmp && mv /outputs/voi/snpEff.config.tmp /outputs/voi/snpEff.config


snpEff build -c /outputs/voi/snpEff.config -gff3 -v mygenome -noCheckCds -noCheckProtein

gunzip -k /outputs/variants/evol1.freebayes.filtered.vcf.gz
snpEff -c /outputs/voi/snpEff.config mygenome /outputs/variants/evol1.freebayes.filtered.vcf > /results/evol1.freebayes.filtered.anno.vcf

gunzip -k /outputs/variants/evol2.freebayes.filtered.vcf.gz
snpEff -c /outputs/voi/snpEff.config -Xmx8g mygenome /outputs/variants/evol2.freebayes.filtered.vcf > /results/evol2.freebayes.filtered.anno.vcf