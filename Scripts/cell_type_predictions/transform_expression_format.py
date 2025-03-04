#!/usr/bin/env python3
"""
Script to transform all cell type prediction CSV files into the format required by PredictDb.

The script:
1. Reads Gene_anno.txt to create a mapping from gene names to Ensembl IDs
2. Finds all *_predictions.csv files in the specified directory
3. For each file:
   a. Changes the header format (first column to NAME)
   b. Replaces gene names with Ensembl IDs
   c. Changes delimiter from comma to tab
   d. Writes the transformed data to a new file in the output directory
"""

import os
import sys
import gzip
import glob
import csv
from pathlib import Path

# Paths - update these to match your file locations
GENE_ANNO_PATH = "/users/PAS2598/duarte63/GitHub/vcf2predictdb-utils/Gene_anno.txt"
INPUT_DIR = "/fs/ess/PAS2598/h5/4_CellType_Preds"
OUTPUT_DIR = "/fs/ess/PAS2598/h5/4_CellType_Preds_formatted"

def create_gene_mapping(anno_file):
    """Create a dictionary mapping gene names to Ensembl IDs from Gene_anno.txt"""
    gene_map = {}
    print(f"Reading gene annotations from {anno_file}")
    
    try:
        with open(anno_file, 'r') as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    ensembl_id = parts[1]
                    gene_name = parts[2]
                    gene_map[gene_name] = ensembl_id
    except Exception as e:
        print(f"Error reading gene annotation file: {e}")
        sys.exit(1)
    
    print(f"Created mapping for {len(gene_map)} genes")
    return gene_map

def transform_data(input_file, output_file, gene_map):
    """Transform the data file to the required format"""
    print(f"Processing {input_file}")
    
    try:
        # Find the actual header line that contains gene_name
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Find the header line
        header_index = -1
        for i, line in enumerate(lines):
            if line.strip().startswith('gene_name'):
                header_index = i
                break
        
        if header_index == -1:
            print(f"Error: Could not find header line starting with 'gene_name' in {input_file}")
            return False
        
        # Extract the data (header and rows)
        header = lines[header_index].strip().split(',')
        data_rows = [line.strip().split(',') for line in lines[header_index+1:] if line.strip()]
        
        # Replace header first column with "NAME"
        header[0] = "NAME"
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Open output file for writing (gzipped)
        with gzip.open(output_file, 'wt') as f_out:
            # Write header
            f_out.write('\t'.join(header) + '\n')
            
            # Process each data row
            missing_genes = []
            for row in data_rows:
                gene_name = row[0]
                # Replace gene name with Ensembl ID if available
                if gene_name in gene_map:
                    row[0] = gene_map[gene_name]
                else:
                    missing_genes.append(gene_name)
                
                # Write the row with tab delimiter
                f_out.write('\t'.join(row) + '\n')
            
            # Report missing genes
            if missing_genes:
                print(f"Warning: Could not find Ensembl IDs for {len(set(missing_genes))} genes in {input_file}:")
                if len(set(missing_genes)) < 10:
                    for gene in set(missing_genes):
                        print(f"  - {gene}")
                else:
                    print(f"  - {', '.join(list(set(missing_genes))[:5])} ... and {len(set(missing_genes))-5} more")
    
    except Exception as e:
        print(f"Error transforming data in {input_file}: {e}")
        return False
    
    print(f"Successfully transformed {input_file} to {output_file}")
    return True

def output_exists(output_file):
    """Check if the output file already exists"""
    return os.path.exists(output_file)

def main():
    # Create gene name to Ensembl ID mapping
    gene_map = create_gene_mapping(GENE_ANNO_PATH)
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find all prediction CSV files
    prediction_files = glob.glob(os.path.join(INPUT_DIR, "*_predictions.csv"))
    
    if not prediction_files:
        print(f"No *_predictions.csv files found in {INPUT_DIR}")
        sys.exit(1)
    
    print(f"Found {len(prediction_files)} prediction files to process")
    
    # Track success/failure/skip
    success_count = 0
    failure_count = 0
    skip_count = 0
    
    # Process each file
    for input_file in prediction_files:
        # Get the base filename and create output path
        base_name = os.path.basename(input_file).replace("_predictions.csv", "")
        output_file = os.path.join(OUTPUT_DIR, f"{base_name}_predictions.txt.gz")
        
        print(f"\nProcessing cell type: {base_name}")
        
        # Skip if output already exists
        if output_exists(output_file):
            print(f"Skipping {base_name} - output file already exists")
            skip_count += 1
            continue
            
        if transform_data(input_file, output_file, gene_map):
            success_count += 1
        else:
            failure_count += 1
    
    print(f"\nSummary: Found {len(prediction_files)} files")
    print(f"  - Success: {success_count}")
    print(f"  - Failure: {failure_count}")
    print(f"  - Skipped: {skip_count}")
    
    if success_count > 0:
        print(f"The transformed files are ready for PredictDb in {OUTPUT_DIR}")

if __name__ == "__main__":
    main()