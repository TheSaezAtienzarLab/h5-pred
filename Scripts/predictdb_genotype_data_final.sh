#!/bin/bash
#SBATCH --job-name=reformat_varids
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=PAS2598
#SBATCH --output=reformat_varids_%j.log
#SBATCH --mail-type=END,FAIL

# Set the directory containing the genotype files
GENOTYPE_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/genotype_files"
# Create a backup directory
BACKUP_DIR="${GENOTYPE_DIR}/original_files_backup"
mkdir -p ${BACKUP_DIR}

echo "Starting variant ID reformatting process at $(date)"
echo "Processing files in: ${GENOTYPE_DIR}"

# Process each genotype file
for GENOTYPE_FILE in ${GENOTYPE_DIR}/chr*_genotypes_final.txt; do
    FILENAME=$(basename "${GENOTYPE_FILE}")
    
    echo "Processing file: ${FILENAME}"
    
    # Create a backup of the original file
    cp "${GENOTYPE_FILE}" "${BACKUP_DIR}/${FILENAME}.bak"
    
    # Use awk to modify the varID column:
    # 1. Extract the header (first line)
    # 2. Process remaining lines to modify varID column
    # 3. Write output to a temporary file
    awk 'NR==1 {print $0; next} 
    {
        # Get the first field (varID)
        varID = $1
        
        # Remove "chr" prefix and add "_GRCh38" suffix
        if (substr(varID, 1, 3) == "chr") {
            newVarID = substr(varID, 4) "_GRCh38"
        } else {
            newVarID = varID "_GRCh38"
        }
        
        # Replace first field and keep others unchanged
        $1 = newVarID
        print $0
    }' "${GENOTYPE_FILE}" > "${GENOTYPE_FILE}.tmp"
    
    # Replace the original file with the modified one
    mv "${GENOTYPE_FILE}.tmp" "${GENOTYPE_FILE}"
    
    echo "Completed processing ${FILENAME}"
done

echo "Reformatting process complete at $(date)"
echo "Original files backed up in: ${BACKUP_DIR}"

# Now also process the variant annotation files to keep them consistent
echo "Processing variant annotation files..."

for ANNOTATION_FILE in ${GENOTYPE_DIR}/chr*_variant_annotation.txt; do
    FILENAME=$(basename "${ANNOTATION_FILE}")
    
    echo "Processing file: ${FILENAME}"
    
    # Create a backup of the original file
    cp "${ANNOTATION_FILE}" "${BACKUP_DIR}/${FILENAME}.bak"
    
    # Modify the variant IDs in the same way
    awk '{
        # Get the first field (varID)
        varID = $1
        
        # Remove "chr" prefix and add "_GRCh38" suffix
        if (substr(varID, 1, 3) == "chr") {
            newVarID = substr(varID, 4) "_GRCh38"
        } else {
            newVarID = varID "_GRCh38"
        }
        
        # Replace first field and keep others unchanged
        $1 = newVarID
        print $0
    }' "${ANNOTATION_FILE}" > "${ANNOTATION_FILE}.tmp"
    
    # Replace the original file with the modified one
    mv "${ANNOTATION_FILE}.tmp" "${ANNOTATION_FILE}"
    
    echo "Completed processing ${FILENAME}"
done

echo "All files processed successfully."