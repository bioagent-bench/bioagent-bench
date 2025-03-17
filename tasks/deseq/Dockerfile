FROM condaforge/mambaforge:24.9.2-0
WORKDIR /app

COPY run_script.sh /app/
COPY run_deseq.R /app/

RUN chmod +x /app/run_script.sh

RUN mamba init

# install some default packages that will be used everywhere
RUN mamba install -y -c bioconda sra-tools r-base

SHELL ["/bin/bash", "-c"]

ENTRYPOINT []

CMD ["/app/run_script.sh"]
