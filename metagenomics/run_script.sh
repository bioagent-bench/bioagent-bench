mkdir data
wget -O data/metagenomics.zip https://zenodo.org/records/7010950/files/dc_workshop.zip
wget -O data/MGRAST_MetaData_JP.xlsx https://zenodo.org/records/7010950/files/MGRAST_MetaData_JP.xlsx

unzip data/metagenomics.zip -d data/

mamba create -n metagenomics
mamba install -c bioconda fastqc