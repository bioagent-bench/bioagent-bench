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
