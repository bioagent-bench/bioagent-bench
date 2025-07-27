FROM condaforge/mambaforge:24.9.2-0

WORKDIR /app

COPY run_script.sh /app/
COPY environment.yml /app/

RUN mamba init bash && \
    . ~/.bashrc && \
    # Create the main environment from environment.yml
    mamba env create -f environment.yml && \
    # Make script executable
    chmod +x /app/run_script.sh


SHELL ["/bin/bash", "--login", "-c"]

CMD ["/app/run_script.sh"] 