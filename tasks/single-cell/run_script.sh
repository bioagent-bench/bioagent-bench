mamba create -n single-cell -f environment.yml
mamba activate single-cell

mkdir -p ./data
mkdir -p ./results

python run_analysis.py