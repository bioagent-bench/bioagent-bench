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
```bash
./prefetch /Users/dionizijefa/Documents/entropic-dev/bio-agent-benchmark/data/SRR9136455/SRR9136455
```

```bash
 ./fastq-dump /Users/dionizijefa/Documents/entropic-dev/bio-agent-benchmark/data/SRR9136455/SRR9136455 -split-files -O /Users/dionizijefa/Documents/entropic-dev/bio-agent-benchmark/data/faster-dump
```

## Virus identification
