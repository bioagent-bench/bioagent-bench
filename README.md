![Funded by Next Gen EU](image.png)

# bioagent-bench
Benchmark for evaluating LLM agents in bioinformatics

# Contents

- **src/**
  - **dataset/**: Code for downloading the input and truth files  
- **tasks/**: Individual bioinformatics tasks  
  - **task_name/**  
  - **Dockerfile**: You can reproduce the truth files with this  
  - **run_script.sh**: Script to run the bioinformatics tools  
  - **run_rscript.R**: Sometimes there is an R downstream analysis  
  - **run_pyscript.R**: Sometimes there is a Python downstream analysis

# General instructions
### CLI: dataset downloader
Use the Click-based CLI in `src/dataset.py` to list tasks and download data, reference files, and results.

- **Install uv (if not installed)**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

- **Create a virtual environment and install dependencies**
```bash
uv sync
```

- **Show help**
```bash
uv run python src/dataset.py --help
```

- **List available tasks**
```bash
uv run python src/dataset.py list-tasks
```

- **Download a single task's data**
```bash
uv run python src/dataset.py download --task giab
```

- **Download multiple tasks**
```bash
uv run python src/dataset.py download --task giab --task deseq --task metagenomics
```

- **Download all tasks**
```bash
uv run python src/dataset.py download --all
```

- **Include reference files**
```bash
uv run python src/dataset.py download --task giab --reference
```

- **Download results only**
```bash
# For specific tasks
uv run python src/dataset.py download --task giab --results --no-data --no-reference

# For all tasks
uv run python src/dataset.py download-all-results
```


Files are downloaded under `tasks/<task_id>/`:
- `data/` for input datasets
- `reference/` for reference data when `--reference` is used
- `results/` for evaluation result files when `--results` or `download-all-results` is used

## Some general ideas about the benchmark
In ideal case you should only care about the src/ folder where you can download the input files, analyze, eval them against "truth" files. Different software, versions, assumptions will produce different results. Where there is consensus actual ground truth like in GIAB (Genome in a bottle) or simulated data we can  produce eval metrics. In other cases we evaluate against correctness of output and try to see if there are at least some overlapping outputs. 
> **DISCLAIMER:** Do not expect "Truth" files or evals to be correct unless explicitly stated.

## What is allowed?
Obviously don't prompt the LLM with the scripts used for reproducing the eval files. We have written two prompts for the humans reading this given below:
1. Data background: Start prompting with this to give some background. This is what it was written
for without reveleaing too much. This is not set in stone, edit it as you wish.
2. Instruct: This is the goal of the analysis. Edit it as you wish.

## Reproducing the pipelines
If you really need to reproduce and check the results to do it you would run
```bash
cd tasks/<task-directory>
```

```bash
docker build -t <task-name> .
```

```bash
docker run -v "$(pwd)/:/app" <task-name>
```

## Datasets and tasks

### Experimental evolution
https://genomics.sschmeier.com/
#### Data background
The experiment follows a similar strategy as in what is called an "experimental evolution" experiment. The final aim is to identify the genome variations in evolved lines of E. coli. The data is composed of a single ancestor line and two evolved lines. The data is from a paired-end sequencing run data from an Illumina HiSeq. This data has been post-processed in two ways already. All sequences that were identified as belonging to the PhiX genome have been removed. Illumina adapters have been removed as well already.
 #### Instruct
 The goal is to find to find and annotate the genome variations in the evolved lines of E.coli. Only output those variants which are shared across both evolved lines


### Cystic fibrosis
https://github.com/jmchilton/SnpEffect/blob/master/papers/Protocol/annotate_example_1.sh
#### Data background
The sample dataset is a simulated dataset for finding the generic cause of Cystic fibrosis. The dataset is real sequencing data from CEPH_1463 dataset provided by the Complete Genomics Diversity Panel. It consists of sequencing of a family: 4 grandparents, 2 parents and 11 siblings. A known Mandelian disease mutation has been added on three siblings, taking care to be consistent with the underlying heplotype structure. The goal is to find the mutation causing the mendalian recessive trait - Cystic Fibrosis. The samples with the observed phenotype are
NA12885, NA12886, NA12879.
Family tree:\
* 1. Siblings: NA12879, NA12880, NA12881, NA12882, NA12883, NA12884, NA12885, NA12886, NA12887 NA12888, NA12893
* 2. Parents: NA12877, NA12878
* 3. Parents of NA12877: NA12889, NA12890
* 4. Parents of NA12878: NA12891, NA12892

### RNA-Seq Differential Expression 
https://docs.tinybio.cloud/docs/rna-seq-tutorial-from-fastq-to-figures
https://pubmed.ncbi.nlm.nih.gov/25233198/
#### Data background
The dataset consists of RNA-Seq samples from Candida parapsilosis wild-type (WT) strains grown in planktonic and biofilm conditions, generated as part of a study on gene expression and biofilm formation. The samples were sequenced on the Illumina HiSeq 2000 platform. The goal of this analysis is to perform differential expression analysis using DESeq2 to identify genes that are significantly up- or down-regulated between planktonic and biofilm conditions, providing insights into biofilm-associated transcriptional changes.

### Transcript Quantification
#### Data background
The dataset contains simulated paired-end RNA-Seq data designed for testing RNA-Seq analysis workflows. The data represents a small set of genes from a human cell line, providing a controlled environment for validating RNA-Seq analysis pipelines. This dataset is derived from the GEUVADIS project, specifically from sample ERR188297, which was processed and simulated to create a manageable test dataset. The test data was generated by selecting 66 expressed genes with read counts exceeding 100,000, resulting in 582 transcripts used to create the transcriptome reference. Transcript-level counts were used to simulate paired-end reads using polyester, with reads shuffled while maintaining pairing and quality scores added using bbmap.

#### Instruct
The goal is to perform transcript quantification from paired-end RNA-Seq reads using the provided transcriptome reference. Analyze the simulated RNA-Seq data (reads_1.fq.gz and reads_2.fq.gz) against the transcriptome reference (transcriptome.fa) to quantify transcript expression levels. The analysis should demonstrate proper RNA-Seq workflow validation using this controlled test dataset containing 582 transcripts from 66 selected genes. The goal is to get the exact matches on the number of reads since the data is simulated.

### Metagenomics
https://jose.theoj.org/papers/10.21105/jose.00209
https://carpentries-lab.github.io/metagenomics-analysis/
#### Data background
The metagenomes that we will use were collected in Cuatro Ciénegas, in a study about the response of the Cuatro Cienegas' bacterial community to nutrient enrichment. In this study, authors compared the differences between the microbial community in its natural, oligotrophic, phosphorus-deficient environment, a pond from the Cuatro Ciénegas Basin (CCB), and the same microbial community under a fertilization treatment. Sample data is Control mecocosm (JC1A) and fertilized pond (JP4D)

### Comparative Genomics
https://www.ahl27.com/CompGenomicsBioc2022/articles/CompGenomicsBioc2022.html

#### Data background
The datasets consists FASTA sequences and GFF annotations of a microbial genome for Micrococcus. The goal of is to do phylogenetic reconstruction of clusters of orthologous co-evolving genes; identify functionally conserved gene clusters across the genomes and group them into co-evolving functional modules. The COGs needs to be filtered based on the following quality criteria. 
1. Clusters present in all 4 organisms
2. Only present in the coding regions
3. Must have at least 1 high confidence annotation
The final result is clustering of the co-evolving genes into functional (annotated clusters) in the form of a .csv file with a "cluster_number" and "consensus_annotation" columns.
The consensus_annotation column may contain a KEGG orthology ID (e.g., K07222) followed by a description.

#### Goal
The final result is clustering of the co-evolving genes into functional (annotated clusters) in the form of a .csv file with a "cluster_number" and "consensus_annotation" columns.
The consensus_annotation column may contain a KEGG orthology ID (e.g., K07222) followed by a description.

Your task is to always use the KEGG orthology ID (the code starting with 'K' followed by numbers) when referencing or summarizing annotations, not the longer description. If there are multiple K numbers, consider all of them but do not invent or merge their meanings beyond what is given.

Example Output:
1,"K07222 K07222, putative flavoprotein involved in K+ transport"

### Viral metagenomcis
https://www.hadriengourle.com/tutorials/metavir/
#### Data background
Real dataset published in 2017 in a study in dolphins, where fecal samples where prepared for viral metagenomics study. The dolphin had a self-limiting gastroenteritis of suspected viral origin. We have available two reads from the fecal samples. We want to find the species of the potential viral agent

#### Goal
I have paired-end sequencing data from a fecal sample of a dolphin that experienced gastroenteritis of suspected viral origin. I need to identify viral species present in this sample. Please analyze the sequencing data and provide the results in a CSV format with the following columns:
contig_count: number of contigs matching each classification
domain: taxonomic domain (e.g., Viruses)
species: name of the viral species

### Alzheimer Mouse Models
Alzheimer's Disease (AD) is a progressive neurodegenerative disorder characterized by hallmark pathologies such as amyloid-β plaques, tau protein tangles, synaptic dysfunction, and neuroinflammation. While mice do not naturally develop Alzheimer's Disease, transgenic mouse models have been developed to recapitulate key pathological features of the disease, enabling researchers to study its underlying mechanisms and evaluate potential therapeutic interventions.

#### Source
This dataset and documentation is sourced from [@evanpeikon/Mouse_Alz_Models](https://github.com/evanpeikon/Mouse_Alz_Models)

#### Data Sources
The data comes from three GEO studies examining different transgenic mouse models of Alzheimer's disease:
- GSE168137: 5xFAD mouse model characterization
- GSE161904: 3xTg-AD mouse model transcriptional profiling
- GSE118523: P301S tau transgenic mice CNS profiling

#### Dataset Files
1. GSE168137_countList.txt (5xFAD)
- Control samples: BL6 wild-type mice (10 total)
  - 5 females: BL6_cortex_8mon_Female_(322,338,340,348,351)
  - 5 males: BL6_cortex_8mon_Male_(389,396,399,410,412)
- Treatment samples: 5xFAD mice (10 total)
  - 5 females: 5xFAD;BL6_cortex_8mon_Female_(295,312,339,341,342)
  - 5 males: 5xFAD;BL6_cortex_8mon_Male_(299,300,307,387,390)
- Age: 8 months
- Tissue: Cortex
- Format: Tab-separated gene counts matrix

2. GSE161904_Raw_gene_counts_cortex.txt (3xTG-AD)
- Control samples: G3R10_Cortex_WT, G3R7_Cortex_WT, G3R9_Cortex_WT
- Treatment samples: G3R1_Cortex_3xTgAD, G3R3_Cortex_3xTgAD, G3R4_Cortex_3xTgAD
- Format: Tab-separated gene counts matrix with ENSEMBL IDs

3.  DEA_PS3O1S.csv (PS3O1S)
- Control samples: Wild-type (WT) mice, 3 replicates (FPKM.tau_01 to tau_03)
- Treatment samples: PS301S transgenic mice, 3 replicates (FPKM.tau_04 to tau_06)
- Format: CSV with differential expression analysis results including:
  - Fold change (fc)
  - Log2 fold change (log2fc)
  - P-value (pval)
  - Q-value (qval)
  - Mean FPKM values for both conditions

#### Goal
Perform a comparative analysis of three different Alzheimer's Disease mouse models (5xFAD, 3xTG-AD, and PS3O1S) to identify shared molecular pathways. The analysis should:

1. Process and normalize the raw gene expression data for each model
2. Perform differential expression analysis with the following criteria:
   - Adjusted p-value < 0.075
   - Absolute log2 fold change > 0.75
3. Conduct pathway enrichment analysis using KEGG pathways for each model
4. Generate a final comparison CSV file that contains:
   - Common pathways present across all three models
   - Corresponding p-values for each pathway-model combination
   - Pathways sorted by minimum p-value across models

The output should be a CSV file with columns: Pathway, 5xFAD_pvalue, 3xTG_AD_pvalue, PS3O1S_pvalue, where each row represents a pathway found in all three models.

### Single Cell RNA Seq
https://github.com/evanpeikon/scRNA-seq-guide
#### Data background
The dataset consists of single-cell RNA sequencing data from human skeletal muscle samples collected before and after acute exercise. The data comes from a study examining cellular responses to exercise in human skeletal muscle, with samples from three subjects collected at pre-exercise and post-exercise timepoints. The sequencing was performed using 10X Genomics technology, generating sparse matrix files (matrix.mtx), cell barcodes (barcodes.tsv), and gene features (features.tsv) for each sample. This dataset enables investigation of how different cell types within skeletal muscle tissue respond to exercise at the single-cell level, revealing cell type-specific transcriptional changes that would be masked in bulk RNA-seq analysis.

#### Goal
The goal is to perform comprehensive single-cell RNA-seq analysis to identify cell types present in human skeletal muscle and characterize their specific responses to acute exercise. The analysis should include quality control, normalization, dimensionality reduction, clustering, cell type identification using marker genes, and differential expression analysis between pre- and post-exercise conditions within each cell type. The final output should be a CSV file containing all differentially expressed genes across cell types, with annotations indicating the predicted cell type for each cluster and the direction of expression changes in response to exercise.

### GIAB Variant Calling
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

#### Data processing
To minimize computational load while preserving mapping fidelity, we first aligned the entire NA12878 HiSeq 300 × read set to the full GRCh38 “no-alt” reference with BWA-MEM. Whole-genome alignment ensures that reads originating from repetitive or paralogous loci retain their correct genomic context, preventing misplacement onto homologous segments of chromosome 20. After alignment we extracted only those read pairs whose primary mapping coordinates fell within chr20 using samtools view, producing a chr20-specific BAM/FASTQ subset (~5 GB) suitable for rapid variant calling. This strategy retains the accuracy of a full-genome workflow while reducing data volume and runtime by more than an order of magnitude, enabling efficient benchmarking against GIAB’s high-confidence chr20 truth set.