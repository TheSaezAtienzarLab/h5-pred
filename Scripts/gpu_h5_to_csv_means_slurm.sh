#SBATCH --account=PAS2598
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28    
#SBATCH --gpus-per-node=2
#SBATCH --mem=32GB
#SBATCH --job-name=h5_test
#SBATCH --output=h5test%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=duarte63@osu.edu

# Load latest miniconda version
module load miniconda3/24.1.2-py310
module load cuda/11.8.0

# Define conda environment path explicitly
export CUDA_VISIBLE_DEVICES=0,1
export CONDA_ENV_PATH="/users/PAS2598/duarte63/.conda/envs/gene_pred"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export PYTHONPATH="$CONDA_ENV_PATH/lib/python3.10/site-packages:$PYTHONPATH"

# Activate conda environment
source activate gene_pred

# Set directories
H5_DIR="/fs/ess/PAS2598/h5/Enformer_output_4bins"
OUTPUT_DIR="/fs/ess/PAS2598/h5/processed_predictions"

# Run the processing script
python gpu_h5_to_csv_means.py $H5_DIR $OUTPUT_DIR