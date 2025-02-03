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

## Virus identification
1. Download the soil paired sequences
```bash
prefetch.sh
```

2. Run the preprocessing scripts
```bash
fastdump.sh
```

3. (Optional) Calculate checksums for the sequences - This is used to create an already supplied input file
```bash
md5sum.sh
```



