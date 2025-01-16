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

```bash
 ./scan/scripts/run_scan10.sh
```

## Virus identification
