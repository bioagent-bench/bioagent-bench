FROM condaforge/mambaforge:latest

WORKDIR /pipeline

COPY /src /pipeline/src
RUN conda env create -f /pipeline/src/environment.yaml -n virus

COPY . /pipeline

SHELL ["/bin/bash", "-c"]

RUN echo "conda activate virus" >> ~/.bashrc

CMD ["sleep", "infinity"] 