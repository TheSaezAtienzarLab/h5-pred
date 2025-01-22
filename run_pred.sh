#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=64GB
#SBATCH --job-name=h5_process
#SBATCH --output=h5_process_%j.log

# Load latest miniconda version
module load miniconda3/24.1.2-py310

# Define conda environment path explicitly
export CONDA_ENV_PATH="/users/PAS2598/duarte63/.conda/envs/gene_pred"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export PYTHONPATH="$CONDA_ENV_PATH/lib/python3.10/site-packages:$PYTHONPATH"

# Activate conda environment
source activate gene_pred

# Run the Python script
python process_h5_files.py