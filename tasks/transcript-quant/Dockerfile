# Use a base image with conda
FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
WORKDIR /app

RUN mkdir -p /app/data
RUN wget -O /app/data/data.tar.gz "https://osf.io/download/684fefc97139f9f5536df7a8/"

RUN cd /app/data && tar -xzf data.tar.gz && rm data.tar.gz

# Copy the script
COPY run_script.sh /app/
COPY environment.yml /app/
RUN chmod +x /app/run_script.sh

# Create and activate the conda environment with salmon
RUN conda env create -f environment.yml && \
    echo "source activate t-quant" > ~/.bashrc

# Make sure conda environment is activated
SHELL ["conda", "run", "-n", "t-quant", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "-n", "t-quant"]
CMD ["./run_script.sh"]