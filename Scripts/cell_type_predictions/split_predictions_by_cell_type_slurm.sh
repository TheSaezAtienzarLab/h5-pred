#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --job-name=reorg_preds
#SBATCH --time=1:00:00            # 6 hours should be sufficient for file operations
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1       # Single task since we're mainly doing I/O
#SBATCH --cpus-per-task=8         # 8 cores should handle parallel file operations
#SBATCH --gpus-per-node=1         # One GPU is enough for this task
#SBATCH --mem=64G                 # 64GB should handle chunked processing
#SBATCH --partition=gpu           
#SBATCH --output=reorg_job_%j.log
#SBATCH --error=reorg_job_%j.err

module load miniconda3/24.1.2-py310
module load cuda/12.4.1

export CUDA_VISIBLE_DEVICES=0,1,2,3
export CONDA_ENV_PATH="/users/PAS2598/duarte63/.conda/envs/gene_pred"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export PYTHONPATH="$CONDA_ENV_PATH/lib/python3.10/site-packages:$PYTHONPATH"
export OMP_NUM_THREADS=24         # Match cpus-per-task
export MKL_NUM_THREADS=24         # For Intel MKL operations

# Run the Python script
source activate gene_pred
python split_predictions_by_cell_type.py