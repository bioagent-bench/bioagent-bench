# Use mambaforge as base image
FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy your scripts and environment files
COPY tasks/deseq/run_script.sh /app/
COPY tasks/deseq/run_deseq.R /app/
COPY tasks/deseq/environment.yml /app/
COPY tasks/deseq/r-environment.yml /app/

# Initialize mamba for shell usage and create environments from yml files
RUN mamba init bash && \
    . ~/.bashrc && \
    # Create the main environment from environment.yml
    mamba env create -f environment.yml && \
    # Create R environment from r-environment.yml
    mamba env create -f r-environment.yml && \
    # Make script executable
    chmod +x /app/run_script.sh

# Set bash as the shell and source bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Default command
CMD ["/app/run_script.sh"]
