mamba create -n comparative
mamba activate comparative
mamba install -c conda-forge -c bioconda r-base r-biocmanager bioconductor-decipher bioconductor-synextend

mkdir -p ./data
mkdir -p ./outputs
mkdir -p ./results

wget -O ./data/micrococcus_data.tar.gz "https://osf.io/download/nz3sq/"
tar -xzf ./data/micrococcus_data.tar.gz -C ./data
rm ./data/micrococcus_data.tar.gz