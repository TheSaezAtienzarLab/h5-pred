import pandas as pd
import numpy as np
import cupy as cp
from pathlib import Path
import os
import glob
import logging
from tqdm import tqdm

def ensure_output_directory(output_dir):
    """Create output directory if it doesn't exist."""
    logging.info(f"Creating output directory: {output_dir}")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

def process_batch(features_batch, region_maps_batch):
    """Process a batch of features using GPU acceleration."""
    try:
        # Extract numeric columns only, excluding 'gene_name' and 'chromo'
        numeric_cols = features_batch.select_dtypes(include=[np.number]).columns
        features_data = features_batch[numeric_cols].values
        
        # Move numeric data to GPU
        features_gpu = cp.array(features_data, dtype=cp.float32)
        results = []
        
        for gene_name, regions in region_maps_batch:
            if regions and isinstance(regions, str):
                region_list = regions.split('|')
                # Get indices of regions in features batch
                region_indices = features_batch['gene_name'].isin(region_list)
                
                if region_indices.any():
                    # Get features for matched regions on GPU
                    region_features = features_gpu[region_indices]
                    # Calculate mean on GPU
                    mean_features = cp.mean(region_features, axis=0)
                    # Move result back to CPU
                    mean_features_cpu = cp.asnumpy(mean_features)
                    results.append([gene_name] + mean_features_cpu.tolist())
        
        return results
    except Exception as e:
        logging.error(f"Error in batch processing: {str(e)}")
        return []

def process_single_file(gene_windows_df, features_path, output_path):
    """Process a single features file using GPU acceleration."""
    # Define numeric column names we expect (0 through 5312)
    numeric_cols = [str(i) for i in range(5313)]
    try:
        logging.info(f"Processing file: {features_path}")
        
        # Skip if output already exists
        if Path(output_path).exists():
            logging.info(f"Skipping {features_path} - output already exists")
            return True
        
        # Read features file in chunks for memory efficiency
        chunk_size = 5000  # Adjust based on your GPU memory
        features_chunks = pd.read_csv(
            features_path,
            chunksize=chunk_size,
            dtype={'gene_name': str, 'chromo': str, **{str(i): np.float32 for i in range(5313)}}
        )
        
        all_results = []
        
        # Process gene windows in batches
        gene_batch_size = 1000
        gene_batches = [
            gene_windows_df[i:i + gene_batch_size] 
            for i in range(0, len(gene_windows_df), gene_batch_size)
        ]
        
        for features_chunk in tqdm(features_chunks, desc="Processing features chunks"):
            for gene_batch in gene_batches:
                batch_results = process_batch(
                    features_chunk,
                    gene_batch[['gene_name', 'region_names']].values
                )
                all_results.extend(batch_results)
        
        # Create final DataFrame
        if all_results:
            result_df = pd.DataFrame(
                all_results,
                columns=['gene_name'] + [str(i) for i in range(5313)]
            )
            
            # Save results using compression
            result_df.to_csv(output_path, index=False, compression='gzip')
            logging.info(f"Successfully processed {len(result_df)} genes for {features_path}")
            return True
        
        return False
    
    except Exception as e:
        logging.error(f"Failed to process file {features_path}: {str(e)}")
        return False

def main():
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('gene_features_processing.log'),
            logging.StreamHandler()
        ]
    )
    
    logging.info("Starting GPU-accelerated gene features processing")
    
    # Define paths
    gene_windows_path = "/users/PAS2598/duarte63/GitHub/h5-pred/gene_window_regions_with_names.csv"
    features_dir = "/fs/ess/PAS2598/h5/processed_predictions"
    output_dir = "/fs/ess/PAS2598/h5/gene_mean_features"
    
    # Ensure output directory exists
    ensure_output_directory(output_dir)
    
    # Read gene windows file
    logging.info(f"Reading gene windows file: {gene_windows_path}")
    gene_windows_df = pd.read_csv(gene_windows_path)
    
    # Pre-filter gene windows to only include valid entries
    gene_windows_df = gene_windows_df[gene_windows_df['region_names'].notna()]
    
    # Get list of all processed feature files
    feature_files = glob.glob(os.path.join(features_dir, "*_processed.csv"))
    logging.info(f"Found {len(feature_files)} feature files to process")
    
    # Process each file
    successful = 0
    failed = 0
    
    for feature_file in tqdm(feature_files, desc="Processing files"):
        output_path = os.path.join(output_dir, Path(feature_file).stem + "_mean_features.csv.gz")
        
        if process_single_file(gene_windows_df, feature_file, output_path):
            successful += 1
        else:
            failed += 1
            
    logging.info("\nProcessing Summary:")
    logging.info(f"Successfully processed: {successful} files")
    logging.info(f"Failed to process: {failed} files")

if __name__ == "__main__":
    main()