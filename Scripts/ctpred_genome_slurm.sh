#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --time=1:00:00           # Reduced to 1 hour for testing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4      # Reduced to 4 cores
#SBATCH --gpus-per-node=1        # Keeping 1 GPU since you likely need it for testing
#SBATCH --mem=32GB               # Reduced to 32GB
#SBATCH --job-name=ctpred_test   # Changed name to indicate test run
#SBATCH --output=predict_test_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=duarte63@osu.edu

module load miniconda3/24.1.2-py310
module load cuda/12.4.1

export CUDA_VISIBLE_DEVICES=0,1,2,3
export CONDA_ENV_PATH="/users/PAS2598/duarte63/.conda/envs/gene_pred"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export PYTHONPATH="$CONDA_ENV_PATH/lib/python3.10/site-packages:$PYTHONPATH"

source activate gene_pred
python ctpred_genome.py