![Funded by Next Gen EU](image.png)

# bio-agent-benchmark

Benchmark for evaluating LLM agents in bioinformatics

## Setup

1. Create a venv with the uv package manager
    ```bash
    uv venv
    ```

2. Init the project inside the venv
    ```bash
    uv init
    ```

### You can install packages using
    ```bash
    uv add jupyterlab
    ```

## Running Bio Workflows

## Experimental evolution
### Data background
 The experiment follows a similar strategy as in what is called an “experimental evolution” experiment. The final aim is to identify the genome variations in evolved lines of E. coli. The data is composed of a single ancestor line and two evolved lines. The data is from a paired-end sequencing run data from an Illumina HiSeq. This data has been post-processed in two ways already. All sequences that were identified as belonging to the PhiX genome have been removed. Illumina adapters have been removed as well already.

## Cystic fibrosis
The sample dataset is a simulated dataset for finding the generic cause of Cystic fibrosis. The dataset is real sequencing data from CEPH_1463 dataset provided by the Complete Genomics Diversity Panel. It consists of sequencing of a family: 4 grandparents, 2 parents and 11 siblings. A known Mandelian disease mutation has been added on three siblings, taking care to be consistent with the underlying heplotype structure. The goal is to find the mutation causing the mendalian recessive trait - Cystic Fibrosis. The samples with the observed phenotype are
NA12885, NA12886, NA12879.
Family tree:\
* 1. Siblings: NA12879, NA12880, NA12881, NA12882, NA12883, NA12884, NA12885, NA12886, NA12887 NA12888, NA12893
* 2. Parents: NA12877, NA12878
* 3. Parents of NA12877: NA12889, NA12890
* 4. Parents of NA12878: NA12891, NA12892
