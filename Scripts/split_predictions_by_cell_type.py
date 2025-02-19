import cudf
import cupy as cp
import os
from pathlib import Path
import glob
import numpy as np

def get_cell_type_columns():
    """Return dictionary mapping cell types to their full column names."""
    return {
        'pericytes': 'pericytes_mean_expression_predicted_expression',
        'schwann': 'schwann_mean_expression_predicted_expression', 
        'lymphocytes': 'lymphocytes_mean_expression_predicted_expression',
        'microglia': 'microglia_mean_expression_predicted_expression',
        'meninges': 'meninges_mean_expression_predicted_expression',
        'astrocytes': 'astrocytes_mean_expression_predicted_expression',
        'ependymal': 'ependymal_mean_expression_predicted_expression',
        'opc': 'opc_mean_expression_predicted_expression',
        'oligodendrocytes': 'oligodendrocytes_mean_expression_predicted_expression', 
        'neurons': 'neurons_mean_expression_predicted_expression',
        'endothelial': 'endothelial_mean_expression_predicted_expression'
    }

def process_predictions_gpu(input_dir, output_dir):
    """Process all prediction files using GPU acceleration."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    cell_type_cols = get_cell_type_columns()
    prediction_files = glob.glob(os.path.join(input_dir, "*_predicted_expressions.csv.gz"))
    print(f"Found {len(prediction_files)} prediction files to process")
    
    # Process files for each cell type
    for cell_type in cell_type_cols.keys():
        print(f"Processing cell type: {cell_type}")
        output_file = os.path.join(output_dir, f"{cell_type}_predictions.csv")
        
        # Create empty output file with headers
        with open(output_file, 'w') as f:
            f.write('sample,' + cell_type_cols[cell_type] + '\n')
        
        # Process each input file with GPU acceleration
        for file_path in prediction_files:
            sample_name = Path(file_path).stem.replace('_predicted_expressions.csv', '')
            print(f"Processing {sample_name} for {cell_type}")
            
            try:
                # Read data directly to GPU memory using cuDF
                needed_cols = [cell_type_cols[cell_type]]
                df = cudf.read_csv(file_path,
                                 compression='gzip',
                                 usecols=needed_cols)
                
                # Add sample column and select needed columns
                df.insert(0, 'sample', sample_name)
                df = df[['sample'] + needed_cols]
                
                # Write results back to disk
                df.to_csv(output_file,
                         mode='a',
                         header=False,
                         index=False)
                
                # Explicitly free GPU memory
                del df
                cp.get_default_memory_pool().free_all_blocks()
                
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                continue

def main():
    # Define directories 
    input_dir = "/fs/ess/PAS2598/h5/predicted_expressions"
    output_dir = "/fs/ess/PAS2598/h5/cell_type_predictions"
    
    # Process all predictions with GPU acceleration
    process_predictions_gpu(input_dir, output_dir)
    print("Done!")

if __name__ == "__main__":
    main()