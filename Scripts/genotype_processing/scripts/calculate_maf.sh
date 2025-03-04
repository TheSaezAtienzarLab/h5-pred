#!/bin/bash
#SBATCH --job-name=maf_calc
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --account=PAS2598
#SBATCH --output=maf_calc_%j.log

# Set paths with updated locations
SCRIPT_DIR="$SLURM_SUBMIT_DIR"
GENOTYPE_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/genotype_files/genotype_files"
ANNOTATION_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/genotype_files/annotation_files"
OUTPUT_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/maf_annotations"

#Conda
module load miniconda3/24.1.2-py310
source activate bcftools_env

# Create output directory
mkdir -p $OUTPUT_DIR

# Run the Python script
echo "Running MAF calculation script..."
python $SCRIPT_DIR/calculate_maf.py \
    --genotype-dir $GENOTYPE_DIR \
    --annotation-dir $ANNOTATION_DIR \
    --output-dir $OUTPUT_DIR
echo "Job completed at $(date)"

