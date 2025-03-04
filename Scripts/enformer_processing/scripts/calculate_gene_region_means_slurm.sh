#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --gpus-per-node=2
#SBATCH --mem=64GB
#SBATCH --job-name=gene_features
#SBATCH --output=gene_features_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=duarte63@osu.edu

# Load required modules
module load miniconda3/24.1.2-py310
module load cuda/11.8.0

# Define conda environment path explicitly
export CUDA_VISIBLE_DEVICES=0
export CONDA_ENV_PATH="/users/PAS2598/duarte63/.conda/envs/gene_pred"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export PYTHONPATH="$CONDA_ENV_PATH/lib/python3.10/site-packages:$PYTHONPATH"

# Activate conda environment
source activate gene_pred

# Run the processing script
python calculate_gene_region_means.py