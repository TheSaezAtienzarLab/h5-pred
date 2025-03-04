#!/usr/bin/env python3
"""
Simple script to calculate MAF for all SNPs and assign a constant R2 value of 1.
Creates annotation files in the required format:
chromosome pos varID ref_vcf alt_vcf R2 MAF rsid rsid_dbSNP150
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Calculate MAF and create annotation files')
    parser.add_argument('--genotype-dir', required=True, 
                        help='Directory containing genotype files')
    parser.add_argument('--annotation-dir', required=True,
                        help='Directory containing annotation files')
    parser.add_argument('--output-dir', required=True,
                        help='Directory to store new annotation files')
    parser.add_argument('--chromosomes', default='all',
                        help='Comma-separated list of chromosomes to process (default: all)')
    
    return parser.parse_args()

def get_chromosomes(genotype_dir, chrom_str):
    """Get list of chromosomes to process"""
    if chrom_str.lower() == 'all':
        # Get all chromosomes from the directory
        chroms = []
        for filename in os.listdir(genotype_dir):
            if filename.startswith('chr') and filename.endswith('_genotypes_final.txt'):
                chrom = filename.split('_')[0][3:]  # Extract '1' from 'chr1_genotypes_final.txt'
                chroms.append(chrom)
        return sorted(chroms, key=lambda x: int(x) if x.isdigit() else x)
    else:
        # Parse comma-separated list
        return chrom_str.split(',')

def process_chromosome(chrom, genotype_dir, annotation_dir, output_dir):
    """
    Process one chromosome:
    1. Read the genotype file
    2. Calculate MAF for each variant
    3. Create new annotation file with MAF and R2=1
    """
    print(f"Processing chromosome {chrom}...")
    
    # Define file paths
    genotype_file = os.path.join(genotype_dir, f"chr{chrom}_genotypes_final.txt")
    annotation_file = os.path.join(annotation_dir, f"chr{chrom}_variant_annotation.txt")
    output_file = os.path.join(output_dir, f"chr{chrom}_annotation.txt")
    
    # Check if input files exist
    if not os.path.exists(genotype_file):
        print(f"Genotype file not found: {genotype_file}")
        return
        
    if not os.path.exists(annotation_file):
        print(f"Annotation file not found: {annotation_file}")
        return
    
    # Read genotype file to calculate MAF
    print(f"Reading genotype file: {genotype_file}")
    
    # Read the file and parse it
    variant_mafs = {}
    
    with open(genotype_file, 'r') as f:
        # Read header to get sample IDs
        header = f.readline().strip().split()
        num_samples = len(header) - 1  # Subtract 1 for the varID column
        
        # Process each variant
        for line_num, line in enumerate(tqdm(f, desc="Calculating MAF")):
            fields = line.strip().split()
            if len(fields) < 2:  # Skip empty lines
                continue
                
            variant_id = fields[0]
            genotypes = fields[1:]
            
            # Count alleles
            total_alleles = 0
            alt_alleles = 0
            
            for gt in genotypes:
                if gt == '0':
                    total_alleles += 2
                elif gt == '1':
                    total_alleles += 2
                    alt_alleles += 1
                elif gt == '2':
                    total_alleles += 2
                    alt_alleles += 2
                # Skip missing values
            
            # Calculate MAF
            if total_alleles > 0:
                maf = alt_alleles / total_alleles
                # Minor allele frequency should be <= 0.5
                if maf > 0.5:
                    maf = 1 - maf
                variant_mafs[variant_id] = round(maf, 5)
            else:
                variant_mafs[variant_id] = np.nan
    
    print(f"Calculated MAF for {len(variant_mafs)} variants")
    
    # Create output file with the required format
    os.makedirs(output_dir, exist_ok=True)
    
    with open(annotation_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Write header
        f_out.write("chromosome\tpos\tvarID\tref_vcf\talt_vcf\tR2\tMAF\trsid\trsid_dbSNP150\n")
        
        # Process each variant in the annotation file
        for line in tqdm(f_in, desc="Creating annotation file"):
            fields = line.strip().split()
            if len(fields) < 5:
                continue
                
            var_id, chr_with_prefix, pos, ref, alt = fields
            
            # Clean chromosome string by removing 'chr' prefix
            chromosome = chr_with_prefix.replace('chr', '')
            
            # Get MAF
            maf = variant_mafs.get(var_id, np.nan)
            maf_str = str(maf) if not np.isnan(maf) else "NA"
            
            # Assign constant R2=1
            r2 = "1.0"
            
            # Use placeholder for rsIDs
            rsid = "."
            rsid_dbsnp150 = "."
            
            # Write to output file
            f_out.write(f"{chromosome}\t{pos}\t{var_id}\t{ref}\t{alt}\t{r2}\t{maf_str}\t{rsid}\t{rsid_dbsnp150}\n")
    
    print(f"Completed chromosome {chrom}")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Get chromosomes to process
    chromosomes = get_chromosomes(args.genotype_dir, args.chromosomes)
    
    print(f"Starting MAF calculation process at {pd.Timestamp.now()}")
    print(f"Processing chromosomes: {', '.join(chromosomes)}")
    
    # Process each chromosome
    for chrom in chromosomes:
        process_chromosome(chrom, args.genotype_dir, args.annotation_dir, args.output_dir)
    
    print(f"MAF calculation completed at {pd.Timestamp.now()}")
    print(f"Output files are in: {args.output_dir}")

if __name__ == "__main__":
    main()