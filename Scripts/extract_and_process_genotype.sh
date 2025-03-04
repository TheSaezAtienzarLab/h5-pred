#!/bin/bash
#SBATCH --job-name=vcf_process
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=PAS2598
#SBATCH --output=vcf_process_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=duarte63@osc.edu

# Set paths
VCF_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/vcf"
SAMPLES_FILE="/users/PAS2598/duarte63/GitHub/vcf2predictdb-utils/samples_to_extract.txt"
OUTPUT_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/genotype_files"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Ensure clean environment
module purge
module load bcftools

# Function to check if output files already exist
output_exists() {
    local CHR=$1
    local GENOTYPE_FILE="$OUTPUT_DIR/chr${CHR}_genotypes_final.txt"
    local ANNOTATION_FILE="$OUTPUT_DIR/chr${CHR}_variant_annotation.txt"
    
    if [ -f "$GENOTYPE_FILE" ] && [ -f "$ANNOTATION_FILE" ]; then
        return 0  # Files exist
    else
        return 1  # Files don't exist
    fi
}

# Function to process a single chromosome
process_chromosome() {
    local CHR=$1
    echo "===== Processing chromosome $CHR ====="
    
    # Check if output files already exist
    if output_exists $CHR; then
        echo "Output files for chromosome $CHR already exist. Skipping..."
        return 0
    fi
    
    # Define VCF file path for this chromosome
    local VCF_FILE="${VCF_DIR}/ALL.chr${CHR}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz"
    
    # Check if VCF file exists
    if [ ! -f "$VCF_FILE" ]; then
        echo "VCF file not found: $VCF_FILE"
        return 1
    fi
    
    echo "Using VCF file: $VCF_FILE"
    
    # Step 1: Extract variant IDs and genotypes for specified samples
    echo "Extracting genotypes..."
    if [ -f "$SAMPLES_FILE" ]; then
        bcftools view -S $SAMPLES_FILE $VCF_FILE | \
          bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' > $OUTPUT_DIR/chr${CHR}_raw_genotypes.txt
    else
        echo "Sample file not found: $SAMPLES_FILE"
        return 1
    fi
    
    # Step 2: Convert phased genotypes to dosage values (0-2)
    echo "Converting genotypes to dosage values..."
    awk '{
      printf "%s", $1;  # Print variant ID
      for(i=2; i<=NF; i++) {
        # Count number of alternative alleles (1s) in the genotype
        if($i == "0|0") printf "\t0";
        else if($i == "0|1" || $i == "1|0") printf "\t1";
        else if($i == "1|1") printf "\t2";
        else printf "\t.";  # Missing or other values
      }
      printf "\n";
    }' $OUTPUT_DIR/chr${CHR}_raw_genotypes.txt > $OUTPUT_DIR/chr${CHR}_dosage_values.txt
    
    # Step 3: Create header with variant ID and sample IDs
    echo "Creating header..."
    echo -n "varID" > $OUTPUT_DIR/header_chr${CHR}.txt
    if [ -f "$SAMPLES_FILE" ]; then
        # Extract sample IDs from the samples file
        cat $SAMPLES_FILE | while read sample; do
            echo -n -e "\t$sample" >> $OUTPUT_DIR/header_chr${CHR}.txt
        done
        echo "" >> $OUTPUT_DIR/header_chr${CHR}.txt
    else
        echo "Sample file not found: $SAMPLES_FILE"
        return 1
    fi
    
    # Step 4: Combine header with dosage data
    echo "Combining header with dosage data..."
    cat $OUTPUT_DIR/header_chr${CHR}.txt $OUTPUT_DIR/chr${CHR}_dosage_values.txt > $OUTPUT_DIR/chr${CHR}_genotypes_final.txt
    
    # Step 5: Create variant annotation file
    echo "Creating variant annotation file..."
    bcftools view $VCF_FILE | \
      bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%CHROM\t%POS\t%REF\t%ALT\n' > $OUTPUT_DIR/chr${CHR}_variant_annotation.txt
    
    # Clean up temporary files
    rm $OUTPUT_DIR/header_chr${CHR}.txt $OUTPUT_DIR/chr${CHR}_raw_genotypes.txt $OUTPUT_DIR/chr${CHR}_dosage_values.txt
    
    echo "Processing complete for chromosome $CHR"
    echo "Output files:"
    echo "  - $OUTPUT_DIR/chr${CHR}_genotypes_final.txt"
    echo "  - $OUTPUT_DIR/chr${CHR}_variant_annotation.txt"
    echo ""
}

# Main processing

# Log start time
echo "Job started at $(date)"
echo "Processing chromosomes 1-22 in $VCF_DIR"

# Process chromosomes 1-22
for CHR in {1..22}; do
    process_chromosome $CHR
done

echo "All chromosomes processed"
echo "Job completed at $(date)"