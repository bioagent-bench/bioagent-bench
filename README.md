![Funded by Next Gen EU](image.png)

# bio-agent-benchmark

Benchmark for evaluating LLM agents in bioinformatics

## Running Bio Workflows
Specific environments for bioinformatics pipelines are solved via mamba which you get via the miniforge installer https://github.com/conda-forge/miniforge

### Experimental evolution
https://genomics.sschmeier.com/
#### Data background
 The experiment follows a similar strategy as in what is called an “experimental evolution” experiment. The final aim is to identify the genome variations in evolved lines of E. coli. The data is composed of a single ancestor line and two evolved lines. The data is from a paired-end sequencing run data from an Illumina HiSeq. This data has been post-processed in two ways already. All sequences that were identified as belonging to the PhiX genome have been removed. Illumina adapters have been removed as well already.

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

### Metagenomics
https://jose.theoj.org/papers/10.21105/jose.00209
https://carpentries-lab.github.io/metagenomics-analysis/
#### Data background
The metagenomes that we will use were collected in Cuatro Ciénegas, in a study about the response of the Cuatro Cienegas’ bacterial community to nutrient enrichment. In this study, authors compared the differences between the microbial community in its natural, oligotrophic, phosphorus-deficient environment, a pond from the Cuatro Ciénegas Basin (CCB), and the same microbial community under a fertilization treatment. Sample data is Control mecocosm (JC1A) and fertilized pond (JP4D)

### Comparative Genomics
https://www.ahl27.com/CompGenomicsBioc2022/articles/CompGenomicsBioc2022.html

#### Data background
The datasets consists FASTA sequences and GFF annotations of a microbial genome for Micrococcus. The main goal. The goal of is to do phylogenetic reconstruction of clusters of orthologous co-evolving genes. The COGs needs to be filtered based on the following quality criteria.
1. No paralogs
2. Clusters present in all 4 organisms
3. Only present in the coding regions
4. Must have at least 1 high confidence annotation
The final result is clustering of the co-evolving genes into functional (annotated clusters)