FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy environment file and install dependencies
COPY environment.yml /app/
RUN mamba env create -f environment.yml

RUN apt-get update && apt-get install -y unzip

# Add conda initialization to bash profile and activate environment
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate fibrosis_env" >> ~/.bashrc

# Make RUN commands use the conda environment
SHELL ["conda", "run", "-n", "fibrosis_env", "/bin/bash", "-c"]

# Copy the analysis script and set permissions
COPY docker_run_script.sh /app/
RUN chmod +x /app/docker_run_script.sh && \
    mkdir -p /app/data

# Set entrypoint to run with conda environment
ENTRYPOINT ["/bin/bash", "-l"]
CMD ["/app/docker_run_script.sh"] 