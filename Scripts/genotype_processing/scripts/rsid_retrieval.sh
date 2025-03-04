#!/bin/bash
#SBATCH --job-name=rsid_retrieval
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --account=PAS2598
#SBATCH --output=rsid_retrieval_%j.log

# Load required modules
module load miniconda3/24.1.2-py310

# Activate conda environment
source activate bcftools_env

# Move to the directory where the script is located
cd $SLURM_SUBMIT_DIR

# Run the Python script with all parameters
python rsid_retrieval.py \
  --input_dir /fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/maf_annotations \
  --output_dir /fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/maf_annotations_with_rsids \
  --dbsnp_vcf /fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/00-All.vcf.gz \
  --pattern "chr*_annotation.txt" \
  --log_level INFO

# Deactivate conda environment when done
conda deactivate

echo "Job completed at $(date)"