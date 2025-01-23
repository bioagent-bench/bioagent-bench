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

Requirements
 a) 
 ```bash
 sudo apt-get install parallel
 ```   

1. 10x single cell
Requirements:
 a) docker
 b) java
 c) git

```bash
./scan/scripts/fast_dump.sh
```

```bash
 ./scan/scripts/fastqc.sh
```

Before running the actual scan 10 pipeline. You should edit nextflow config files to reflect your system
located in scan10/src/nextflow.config

```bash
 ./scan/scripts/run_scan10.sh
```

## Virus identification
Download the soil paired sequences
```bash
prefetch SRX10854940 SRX10854939 SRX10854938 SRX10854937 SRX10854936 SRX10854935 SRX10854934 SRX10854933
```
