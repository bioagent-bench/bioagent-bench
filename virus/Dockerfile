FROM condaforge/mambaforge:latest

WORKDIR /pipeline
COPY . /pipeline

RUN mamba env create -f /pipeline/environment.yaml
ENV PATH /opt/conda/envs/$(head -1 environment.yaml | cut -d' ' -f2)/bin:$PATH

ENTRYPOINT ["snakemake"]
CMD ["--help"] 