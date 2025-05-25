mamba create -f environment.yml
mkdir -p ./data

# Download the data 
wget -O ./data/alzheimer_mouse_data.tar.gz "https://osf.io/download/6832e9df21129298ddd6a3f8/"
tar -xzf ./data/alzheimer_mouse_data.tar.gz -C ./data --strip-components=1

# Remove the tar.gz file to save space
rm ./data/alzheimer_mouse_data.tar.gz

python run_analysis.py