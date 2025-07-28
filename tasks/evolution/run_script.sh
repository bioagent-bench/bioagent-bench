#!/bin/bash
set -e

# Set dynamic thread count
THREADS=$(nproc)

# Logging functions
TOTAL_STEPS=12
CURRENT_STEP=0

log_step() {
    CURRENT_STEP=$((CURRENT_STEP + 1))
    echo "=================================================================="
    echo "STEP $CURRENT_STEP/$TOTAL_STEPS: $1"
    echo "=================================================================="
    echo "Started at: $(date)"
    echo ""
}

log_substep() {
    echo "  → $1"
}

log_completion() {
    echo ""
    echo "Completed at: $(date)"
    echo ""
}

env_exists() {
    local env_name=$1
    mamba env list | grep -q "^${env_name} "
}

# STEP 1: Environment Setup
log_step "Setting up conda environments"

# Sometimes these environments get spastic and this is easier
if ! env_exists "fastp_env"; then
    log_substep "Creating fastp environment"
    mamba create -n fastp_env bioconda::fastp -y
fi

if ! env_exists "fastqc_env"; then
    log_substep "Creating fastqc environment"
    mamba create -n fastqc_env bioconda::fastqc -y
fi

if ! env_exists "multiqc_env"; then
    log_substep "Creating multiqc environment"
    mamba create -n multiqc_env bioconda::multiqc -y
fi

if ! env_exists "spades_env"; then
    log_substep "Creating spades environment"
    mamba create -n spades_env bioconda::spades -y
fi

if ! env_exists "quast_env"; then
    log_substep "Creating quast environment"
    mamba create -n quast_env bioconda::quast -y
fi

if ! env_exists "mapping_env"; then
    log_substep "Creating mapping environment"
    mamba create -y -n mapping_env -c bioconda samtools bwa
fi

if ! env_exists "variant_env"; then
    log_substep "Creating variant calling environment"
    mamba create -y -n variant_env -c bioconda samtools bamtools freebayes
fi

if ! env_exists "compleasm_env"; then
    log_substep "Creating compleasm environment"
    mamba create -y -n compleasm_env bioconda::compleasm
fi

if ! env_exists "prokka_env"; then
    log_substep "Creating prokka environment"
    mamba create -y -n prokka_env bioconda::prokka
fi

if ! env_exists "voi_env"; then
    log_substep "Creating variant annotation environment"
    mamba create -y -n voi_env -c bioconda snpeff genometools-genometools
fi

if ! env_exists "bcf_env"; then
    log_substep "Creating bcf tools environment"
    mamba create -y -n bcf_env -c bioconda bcftools
fi

log_completion

# Create directory structure
mkdir -p data
mkdir -p outputs/trimmed
mkdir -p outputs/trimmed-fastqc
mkdir -p outputs/assembly
mkdir -p outputs/mappings
mkdir -p outputs/variants
mkdir -p outputs/annotation
mkdir -p outputs/voi
mkdir -p outputs/bcf
mkdir -p results

# Download the data to data folder
cd data
wget -O data.tar.gz https://osf.io/2jc4a/download
tar -xvzf data.tar.gz --strip-components=1
cd ..

# STEP 2: Quality Control - Adapter Trimming
log_step "Quality control - trimming adapters"

log_substep "Trimming adapters for ancestral sample"
mamba run -n fastp_env fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html ./outputs/trimmed/anc.fastp.html --json ./outputs/trimmed/anc.fastp.json -i ./data/anc_R1.fastq.gz \
 -I ./data/anc_R2.fastq.gz -o ./outputs/trimmed/anc_R1.fastq.gz -O ./outputs/trimmed/anc_R2.fastq.gz

log_substep "Trimming adapters for evolved sample 1"
mamba run -n fastp_env fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html ./outputs/trimmed/evol1.fastp.html --json ./outputs/trimmed/evol1.fastp.json -i ./data/evol1_R1.fastq.gz \
 -I ./data/evol1_R2.fastq.gz -o ./outputs/trimmed/evol1_R1.fastq.gz -O ./outputs/trimmed/evol1_R2.fastq.gz

log_substep "Trimming adapters for evolved sample 2"
mamba run -n fastp_env fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $THREADS \
 --html ./outputs/trimmed/evol2.fastp.html --json ./outputs/trimmed/evol2.fastp.json -i ./data/evol2_R1.fastq.gz \
 -I ./data/evol2_R2.fastq.gz -o ./outputs/trimmed/evol2_R1.fastq.gz -O ./outputs/trimmed/evol2_R2.fastq.gz

log_completion

# STEP 3: Quality Control - FastQC Analysis
log_step "Quality control - FastQC analysis"

log_substep "Running FastQC on trimmed reads"
mamba run -n fastqc_env fastqc -t $THREADS -o ./outputs/trimmed-fastqc ./outputs/trimmed/*.fastq.gz

log_completion

# STEP 4: Quality Control - MultiQC Visualization
log_step "Quality control - MultiQC visualization"

log_substep "Generating MultiQC report"
mamba run -n multiqc_env multiqc ./outputs/trimmed-fastqc ./outputs/trimmed -o ./outputs

log_completion

# STEP 5: Genome Assembly
log_step "Genome assembly with SPAdes"

log_substep "Running SPAdes assembly"
mamba run -n spades_env spades.py -o ./outputs/assembly --careful -t $THREADS -1 ./outputs/trimmed/anc_R1.fastq.gz -2 ./outputs/trimmed/anc_R2.fastq.gz

log_completion

# STEP 6: Assembly Quality Assessment
log_step "Assembly quality assessment with QUAST"

log_substep "Running QUAST quality assessment"
mamba run -n quast_env quast -t $THREADS -o ./outputs/assembly/quast ./outputs/assembly/scaffolds.fasta

log_completion

# STEP 7: Read Mapping
log_step "Read mapping to reference genome"

log_substep "Preparing reference genome"
cp ./outputs/assembly/scaffolds.fasta ./outputs/mappings/scaffolds.fasta
mamba run -n mapping_env bwa index ./outputs/mappings/scaffolds.fasta

# log_substep "Mapping evolved sample 1 reads"
mamba run -n mapping_env bash -c "bwa mem -t $THREADS ./outputs/mappings/scaffolds.fasta ./outputs/trimmed/evol1_R1.fastq.gz ./outputs/trimmed/evol1_R2.fastq.gz > ./outputs/mappings/evol1.sam"

# log_substep "Mapping evolved sample 2 reads"
mamba run -n mapping_env bash -c "bwa mem -t $THREADS ./outputs/mappings/scaffolds.fasta ./outputs/trimmed/evol2_R1.fastq.gz ./outputs/trimmed/evol2_R2.fastq.gz > ./outputs/mappings/evol2.sam"

log_completion

# STEP 8: Post-processing Mappings
log_step "Post-processing mapping files"

log_substep "Processing evolved sample 1 mappings"
# Evol1 line
mamba run -n mapping_env bash -c "samtools sort -@ $THREADS -n -O sam ./outputs/mappings/evol1.sam | samtools fixmate -m -O bam - ./outputs/mappings/evol1.fixmate.bam"
mamba run -n mapping_env samtools sort -@ $THREADS -O bam -o ./outputs/mappings/evol1.sorted.bam ./outputs/mappings/evol1.fixmate.bam
mamba run -n mapping_env samtools markdup -@ $THREADS -r -S ./outputs/mappings/evol1.sorted.bam ./outputs/mappings/evol1.sorted.dedup.bam
mamba run -n mapping_env bash -c "samtools view -@ $THREADS -h -b -q 20 ./outputs/mappings/evol1.sorted.dedup.bam > ./outputs/mappings/evol1.sorted.dedup.q20.bam"
mamba run -n mapping_env bash -c "samtools view -@ $THREADS -b -f 4 ./outputs/mappings/evol1.sorted.dedup.bam > ./outputs/mappings/evol1.sorted.unmapped.bam"
rm ./outputs/mappings/evol1.sam
rm ./outputs/mappings/evol1.fixmate.bam
rm ./outputs/mappings/evol1.sorted.bam
rm ./outputs/mappings/evol1.sorted.dedup.bam

log_substep "Processing evolved sample 2 mappings"
# Evol2 line
mamba run -n mapping_env bash -c "samtools sort -@ $THREADS -n -O sam ./outputs/mappings/evol2.sam | samtools fixmate -m -O bam - ./outputs/mappings/evol2.fixmate.bam"
mamba run -n mapping_env samtools sort -@ $THREADS -O bam -o ./outputs/mappings/evol2.sorted.bam ./outputs/mappings/evol2.fixmate.bam
mamba run -n mapping_env samtools markdup -@ $THREADS -r -S ./outputs/mappings/evol2.sorted.bam ./outputs/mappings/evol2.sorted.dedup.bam
mamba run -n mapping_env bash -c "samtools view -@ $THREADS -h -b -q 20 ./outputs/mappings/evol2.sorted.dedup.bam > ./outputs/mappings/evol2.sorted.dedup.q20.bam"
mamba run -n mapping_env bash -c "samtools view -@ $THREADS -b -f 4 ./outputs/mappings/evol2.sorted.dedup.bam > ./outputs/mappings/evol2.sorted.unmapped.bam"

rm ./outputs/mappings/evol2.sam
rm ./outputs/mappings/evol2.fixmate.bam
rm ./outputs/mappings/evol2.sorted.bam
rm ./outputs/mappings/evol2.sorted.dedup.bam

log_completion

# STEP 9: Variant Calling
log_step "Variant calling with FreeBayes"

log_substep "Preparing reference genome index"
mamba run -n mapping_env samtools faidx ./outputs/assembly/scaffolds.fasta

log_substep "Calling variants for evolved sample 1"
mamba run -n variant_env bamtools index -in ./outputs/mappings/evol1.sorted.dedup.q20.bam
# p 1 for bacteria
mamba run -n variant_env bash -c "freebayes -p 1 -f ./outputs/assembly/scaffolds.fasta ./outputs/mappings/evol1.sorted.dedup.q20.bam > ./outputs/variants/evol1.freebayes.vcf"
mamba run -n variant_env bgzip -f ./outputs/variants/evol1.freebayes.vcf
mamba run -n variant_env tabix -p vcf ./outputs/variants/evol1.freebayes.vcf.gz
mamba run -n variant_env bash -c "zcat ./outputs/variants/evol1.freebayes.vcf.gz | vcffilter -f \"QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1\" | bgzip > ./outputs/variants/evol1.freebayes.filtered.vcf.gz"
mamba run -n variant_env tabix -p vcf ./outputs/variants/evol1.freebayes.filtered.vcf.gz

log_substep "Calling variants for evolved sample 2"
mamba run -n variant_env bamtools index -in ./outputs/mappings/evol2.sorted.dedup.q20.bam
mamba run -n variant_env bash -c "freebayes -p 1 -f ./outputs/assembly/scaffolds.fasta ./outputs/mappings/evol2.sorted.dedup.q20.bam > ./outputs/variants/evol2.freebayes.vcf"
mamba run -n variant_env bgzip -f ./outputs/variants/evol2.freebayes.vcf
mamba run -n variant_env tabix -p vcf ./outputs/variants/evol2.freebayes.vcf.gz
mamba run -n variant_env bash -c "zcat ./outputs/variants/evol2.freebayes.vcf.gz | vcffilter -f \"QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1\" | bgzip > ./outputs/variants/evol2.freebayes.filtered.vcf.gz"
mamba run -n variant_env tabix -p vcf ./outputs/variants/evol2.freebayes.filtered.vcf.gz

log_completion

# STEP 10: Genome Annotation
log_step "Genome annotation"

log_substep "Running BUSCO completeness assessment"
mamba run -n compleasm_env compleasm run -a ./outputs/assembly/scaffolds.fasta -o ./outputs/annotation -l bacteria_odb10 -t $THREADS

log_substep "Running Prokka gene annotation"
mamba run -n prokka_env prokka --cpus $THREADS --kingdom Bacteria --genus Escherichia --species coli --outdir ./outputs/annotation/prokka ./outputs/assembly/scaffolds.fasta

log_substep "Preparing SnpEff annotation database"
mkdir -p ./outputs/voi/data/mygenome
cp ./outputs/assembly/scaffolds.fasta ./outputs/voi/data/mygenome/sequences.fa
gzip ./outputs/voi/data/mygenome/sequences.fa
cp ./outputs/annotation/prokka/*.gff ./outputs/voi/data/mygenome/genes.gff
gzip ./outputs/voi/data/mygenome/genes.gff

SNPEFF_CONFIG=$(mamba run -n voi_env bash -c "find ~ -name snpEff.config | head -n 1")
cp "$SNPEFF_CONFIG" ./outputs/voi/snpEff.config
sed -i 's|^data.dir.*|data.dir = data/|' ./outputs/voi/snpEff.config

# Add custom genome configuration to snpEff.config
cat >> ./outputs/voi/snpEff.config << EOF
# Custom genome configuration
mygenome.genome : EColiMut
mygenome.codonTable : Bacterial_and_Plant_Plastid
mygenome.checkProtein : false
mygenome.checkCds : false
EOF

log_substep "Building SnpEff database"
mamba run -n voi_env snpEff build -c ./outputs/voi/snpEff.config -gff3 -v mygenome -noCheckCds -noCheckProtein

log_completion

# STEP 11: Variant Annotation
log_step "Variant annotation with SnpEff"

log_substep "Annotating variants for evolved sample 1"
gunzip -k ./outputs/variants/evol1.freebayes.filtered.vcf.gz
mamba run -n voi_env bash -c "snpEff -c ./outputs/voi/snpEff.config -Xmx16g mygenome ./outputs/variants/evol1.freebayes.filtered.vcf > ./outputs/variants/evol1.freebayes.filtered.anno.vcf"

log_substep "Annotating variants for evolved sample 2"
gunzip -k ./outputs/variants/evol2.freebayes.filtered.vcf.gz
mamba run -n voi_env bash -c "snpEff -c ./outputs/voi/snpEff.config -Xmx16g mygenome ./outputs/variants/evol2.freebayes.filtered.vcf > ./outputs/variants/evol2.freebayes.filtered.anno.vcf"

log_completion

# STEP 12: Final Analysis and Results Generation
log_step "Generating final results"

log_substep "Comparing variants between evolved samples"
mamba run -n bcf_env bash -c '
  bcftools view -Oz -o ./outputs/variants/evol1.vcf.gz \
                ./outputs/variants/evol1.freebayes.filtered.anno.vcf  && \
  bcftools index -t ./outputs/variants/evol1.vcf.gz
'

mamba run -n bcf_env bash -c '
  bcftools view -Oz -o ./outputs/variants/evol2.vcf.gz \
                ./outputs/variants/evol2.freebayes.filtered.anno.vcf  && \
  bcftools index -t ./outputs/variants/evol2.vcf.gz
'

mamba run -n bcf_env bcftools isec -p ./outputs/bcf -Oz \
  ./outputs/variants/evol1.vcf.gz \
  ./outputs/variants/evol2.vcf.gz

log_substep "Extracting shared high-impact variants"
# Create shared variants CSV with high/moderate impact mutations
mamba run -n bcf_env bash -c '
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n" \
               ./outputs/bcf/0002.vcf.gz \
| awk -F"\t" '\''BEGIN{OFS=","; print "CHROM,POS,REF,ALT,GENE,IMPACT,EFFECT,STATUS"}
  {
    split($5, ann, "|");           #  ann[3] = IMPACT, ann[4] = Gene_Name, ann[2] = Annotation/EFFECT
    impact = ann[3];
    gene   = ann[4];
    effect = ann[2];
    if ((impact=="HIGH" || impact=="MODERATE") && !seen[gene]++) {
        print $1,$2,$3,$4,gene,impact,effect,"shared";
    }
  }
'\'' > ./results/variants_shared.csv
'

log_substep "Generating gene annotations for variants"
# Create detailed gene annotations CSV
echo "Gene_Name,COG,Function" > ./results/gene_annotations.csv
if [ -s ./results/variants_shared.csv ] && [ $(wc -l < ./results/variants_shared.csv) -gt 1 ]; then
    grep -Ff <(tail -n +2 ./results/variants_shared.csv | cut -d, -f5) \
             ./outputs/annotation/prokka/*.tsv | \
     awk -F'\t' 'BEGIN { OFS="," } 
    {
        if (NF < 7) next
        gene_name = $4
        cog = $6
        function_desc = $7
        if (gene_name == "") gene_name = "unknown"
        if (cog == "") cog = "no_COG"
        if (function_desc == "") function_desc = "hypothetical protein"
        gsub(/,/, ";", function_desc)
        gsub(/"/, "", function_desc)
        print gene_name, cog, function_desc
    }' >> ./results/gene_annotations.csv
fi