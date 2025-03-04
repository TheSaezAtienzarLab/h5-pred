#!/usr/bin/env python3
import os
import sys
import subprocess
import pandas as pd
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor
import time
import logging
from datetime import timedelta

def setup_logging(output_dir, log_level=logging.INFO):
    """Set up logging to both console and file."""
    # Create logs directory if it doesn't exist
    logs_dir = os.path.join(output_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Create a unique log filename with timestamp
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    log_file = os.path.join(logs_dir, f"rsid_retrieval_{timestamp}.log")
    
    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(log_level)
    
    # Clear any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', 
                                      datefmt='%H:%M:%S')
    console_handler.setFormatter(console_format)
    
    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)
    
    # Add handlers to logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    logging.info(f"Logging initialized. Log file: {log_file}")
    return logger

def process_chromosome(input_file, output_file, dbsnp_vcf):
    """Process a single chromosome file to retrieve RSIDs."""
    logger = logging.getLogger()
    logger.info(f"Processing {input_file}...")
    start_time = time.time()
    
    # Read the annotation file
    df = pd.read_csv(input_file, sep='\t')
    
    # Initialize new rsid column if it doesn't exist
    if 'rsid' not in df.columns:
        df['rsid'] = '.'
    
    # Get chromosome number from the first row
    chr_num = str(df.iloc[0]['chromosome'])
    # Remove "chr" prefix if present for consistency
    chr_num = chr_num.replace('chr', '')
    
    # Track statistics
    total_variants = len(df)
    rsids_found = 0
    checkpoint_interval = 10000  # Save progress every 10,000 variants
    
    logger.info(f"Chr {chr_num}: Starting processing of {total_variants} variants")
    
    # Process each variant
    for idx, row in enumerate(df.itertuples(index=True), 1):
        # Check if we already have an RSID
        if hasattr(row, 'rsid') and row.rsid != '.':
            continue
            
        # Skip variants that don't have position or allele information
        if not hasattr(row, 'pos') or pd.isna(row.pos) or not hasattr(row, 'ref_vcf') or not hasattr(row, 'alt_vcf'):
            continue
            
        pos = int(row.pos)
        ref = row.ref_vcf
        alt = row.alt_vcf
        
        # Query dbSNP for this position
        cmd = f"tabix {dbsnp_vcf} {chr_num}:{pos}-{pos}"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
            
            # If first query fails, try with chr prefix
            if result.returncode != 0 or not result.stdout.strip():
                cmd = f"tabix {dbsnp_vcf} chr{chr_num}:{pos}-{pos}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
            
            if result.returncode == 0 and result.stdout.strip():
                # Parse dbSNP output
                for line in result.stdout.strip().split('\n'):
                    fields = line.split('\t')
                    if len(fields) < 5:
                        continue
                        
                    dbsnp_chr = fields[0].replace('chr', '')
                    try:
                        dbsnp_pos = int(fields[1])
                        dbsnp_rsid = fields[2]
                        dbsnp_ref = fields[3]
                        dbsnp_alt = fields[4]
                    except (ValueError, IndexError):
                        continue
                    
                    # Match by position and alleles
                    if dbsnp_pos == pos and dbsnp_ref == ref and (dbsnp_alt == alt or alt in dbsnp_alt.split(',')):
                        # Update both rsid and rsid_dbSNP150 columns
                        df.at[row.Index, 'rsid'] = dbsnp_rsid
                        if 'rsid_dbSNP150' in df.columns:
                            df.at[row.Index, 'rsid_dbSNP150'] = dbsnp_rsid
                        rsids_found += 1
                        break
        except Exception as e:
            logger.error(f"Error processing position {pos} on chr{chr_num}: {e}")
        
        # Print progress and save checkpoint periodically
        if idx % checkpoint_interval == 0:
            elapsed = time.time() - start_time
            percent_complete = (idx / total_variants) * 100
            est_total_time = elapsed / (idx / total_variants)
            est_remaining = est_total_time - elapsed
            
            elapsed_str = str(timedelta(seconds=int(elapsed)))
            remaining_str = str(timedelta(seconds=int(est_remaining)))
            
            logger.info(f"Chr {chr_num}: {idx:,}/{total_variants:,} variants ({percent_complete:.2f}%) | "
                       f"Found {rsids_found:,} RSIDs | "
                       f"Elapsed: {elapsed_str} | Remaining: {remaining_str}")
            
            # Save checkpoint
            checkpoint_file = output_file.replace('.txt', f'_checkpoint_{idx}.txt')
            df.iloc[:idx].to_csv(checkpoint_file, sep='\t', index=False)
    
    # Save the final dataframe
    df.to_csv(output_file, sep='\t', index=False)
    
    elapsed = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))
    logger.info(f"Chr {chr_num}: Completed in {elapsed_str} | "
               f"Found {rsids_found:,}/{total_variants:,} RSIDs ({(rsids_found/total_variants)*100:.2f}%) | "
               f"Saved to {output_file}")
    return output_file, rsids_found

def main():
    parser = argparse.ArgumentParser(description='Retrieve RSIDs for variants from dbSNP for all chromosomes')
    parser.add_argument('--input_dir', default='/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/maf_annotations',
                       help='Directory containing annotation files')
    parser.add_argument('--output_dir', default='/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/maf_annotations_with_rsids',
                       help='Directory to save results')
    parser.add_argument('--dbsnp_vcf', default='/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/00-All.vcf.gz',
                       help='Path to dbSNP VCF file')
    parser.add_argument('--threads', type=int, default=4, 
                       help='Number of chromosomes to process in parallel')
    parser.add_argument('--pattern', default='chr*_annotation.txt',
                       help='Pattern to match annotation files')
    parser.add_argument('--log_level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Logging level')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(args.output_dir, log_level)
    
    # Log script start with arguments
    logger.info("=" * 80)
    logger.info("RSID Retrieval Script Started")
    logger.info(f"Arguments: {vars(args)}")
    
    # Get list of chromosome files
    input_pattern = os.path.join(args.input_dir, args.pattern)
    chr_files = glob.glob(input_pattern)
    
    if not chr_files:
        logger.error(f"No chromosome annotation files found matching pattern: {input_pattern}")
        return
    
    logger.info(f"Found {len(chr_files)} chromosome annotation files")
    
    # Process chromosomes in parallel
    start_time = time.time()
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for input_file in chr_files:
            # Create output filename
            base_name = os.path.basename(input_file)
            output_name = base_name.replace('_annotation.txt', '_with_rsids.txt')
            output_file = os.path.join(args.output_dir, output_name)
            
            # Submit the task
            logger.info(f"Submitting job for {base_name}")
            future = executor.submit(process_chromosome, input_file, output_file, args.dbsnp_vcf)
            futures.append((future, base_name))
        
        # Wait for tasks to complete and report results
        completed = 0
        for future, base_name in futures:
            try:
                output_file, rsids_found = future.result()
                completed += 1
                logger.info(f"Completed {base_name} ({completed}/{len(chr_files)}): Found {rsids_found:,} RSIDs")
            except Exception as e:
                logger.error(f"Error processing {base_name}: {e}", exc_info=True)
    
    # Log completion summary
    elapsed = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))
    logger.info(f"All chromosomes processed in {elapsed_str}. Results are in {args.output_dir}")
    logger.info("=" * 80)

if __name__ == "__main__":
    main()