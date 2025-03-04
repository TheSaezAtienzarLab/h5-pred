import h5py
import cupy as cp
import numpy as np
import pandas as pd
from pathlib import Path
import os
import sys
from tqdm import tqdm

def process_h5_file(h5_path, output_dir):
    """
    Process H5 file where keys are genes and values are 4x5313 matrices
    Using GPU acceleration for numerical operations only
    """
    try:
        # Create output filename based on input filename
        output_name = Path(h5_path).stem + "_processed.csv"
        output_path = Path(output_dir) / output_name
        
        # Skip if output file already exists
        if output_path.exists():
            print(f"Skipping {h5_path} - output already exists")
            return
        
        gene_features = {}
        with h5py.File(h5_path, 'r') as f:
            genes = list(f.keys())
            print(f"Processing {h5_path} - Total genes found: {len(genes)}")
            
            # Process genes in batches for better GPU utilization
            batch_size = 1000
            for i in tqdm(range(0, len(genes), batch_size), desc="Processing gene batches"):
                batch_genes = genes[i:i + batch_size]
                batch_data = []
                
                # Collect batch data
                for gene in batch_genes:
                    data = f[gene][:]
                    batch_data.append(data)
                
                if batch_data:
                    # Move data to GPU for processing
                    batch_array = cp.array(batch_data)
                    # Calculate means on GPU
                    mean_features = cp.mean(batch_array, axis=1)
                    # Move results back to CPU
                    mean_features_cpu = cp.asnumpy(mean_features)
                    
                    # Store results using CPU operations
                    for gene, features in zip(batch_genes, mean_features_cpu):
                        gene_features[gene] = features.tolist()
        
        # Create DataFrame on CPU
        df = pd.DataFrame.from_dict(gene_features, orient='index')
        
        # Add gene name column
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'gene_name'}, inplace=True)
        
        # Extract chromosome using pandas string operations
        df['chromo'] = df['gene_name'].str.split('_').str[0]
        
        # Save to CSV
        df.to_csv(output_path, index=False)
        print(f"Saved processed data to {output_path}")
        print(f"Shape: {df.shape}")
        
    except Exception as e:
        print(f"Error processing file {h5_path}: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_h5.py <h5_directory> <output_directory>")
        sys.exit(1)
        
    h5_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of all H5 files
    h5_files = list(Path(h5_dir).glob("*.h5"))
    print(f"Found {len(h5_files)} H5 files to process")
    
    # Process each file
    for h5_path in h5_files:
        process_h5_file(str(h5_path), output_dir)

if __name__ == "__main__":
    main()