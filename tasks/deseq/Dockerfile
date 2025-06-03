FROM condaforge/mambaforge:24.9.2-0

WORKDIR /app

# Copy necessary files
COPY run_script.sh /app/
COPY run_deseq.R /app/
COPY environment.yml /app/
COPY r-environment.yml /app/

# Make the run script executable
RUN chmod +x /app/run_script.sh

# Initialize mamba
RUN mamba init bash

# Create the RNA-seq environment and install dependencies
RUN mamba env create -f environment.yml && \
    mamba env create -f r-environment.yml && \
    mamba clean -a -y

# Install system dependencies
RUN apt-get update && \
    apt-get install -y procps wget curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set up R environment
RUN mamba run -n rnaseq-r R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(version = "3.18", ask = FALSE); BiocManager::install(c("DESeq2", "BiocParallel", "clusterProfiler", "GenomeInfoDbData", "GO.db"), ask = FALSE); if (!require("ggpubr")) install.packages("ggpubr", dependencies = TRUE)'

# Set the shell to bash
SHELL ["conda", "run", "-n", "rnaseq", "/bin/bash", "-c"]

# Set the entrypoint to use bash with login shell (to ensure conda initialization)
ENTRYPOINT ["/bin/bash", "-l", "-c"]

# Default command
# CMD ["/app/run_script.sh"]
