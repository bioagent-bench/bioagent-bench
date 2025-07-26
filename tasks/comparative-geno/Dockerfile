FROM condaforge/mambaforge:24.9.2-0

# Set working directory
WORKDIR /app

# Copy your scripts and environment files
COPY run_script.sh /app/
COPY run_comparative.R /app/
COPY environment.yml /app/

# Initialize mamba for shell usage and create environments from yml files
RUN mamba init bash && \
    . ~/.bashrc && \
    # Create the main environment from environment.yml
    mamba env create -f environment.yml && \
    # Make script executable
    chmod +x /app/run_script.sh

# Set bash as the shell and ensure mamba is available
SHELL ["/bin/bash", "--login", "-c"]


CMD ["mamba", "run", "-n", "comparative", "/app/run_script.sh"]