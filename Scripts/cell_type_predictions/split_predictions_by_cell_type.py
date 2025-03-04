import os
from pathlib import Path
import glob
import numpy as np
import pandas as pd

def get_cell_type_columns():
    """Return dictionary mapping cell types to their full column names."""
    return {
        'pericytes': 'pericytes_predicted',
        'schwann': 'schwann_predicted', 
        'lymphocytes': 'lymphocytes_predicted',
        'microglia': 'microglia_predicted',
        'meninges': 'meninges_predicted',
        'astrocytes': 'astrocytes_predicted',
        'ependymal': 'ependymal_predicted',
        'opc': 'opc_predicted',
        'oligodendrocytes': 'oligodendrocytes_predicted', 
        'neurons': 'neurons_predicted',
        'endothelial': 'endothelial_predicted'
    }

def process_predictions(input_dir, output_dir):
    """Process all prediction files to create gene x sample matrices."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    cell_type_cols = get_cell_type_columns()
    prediction_files = glob.glob(os.path.join(input_dir, "*.csv"))
    print(f"Found {len(prediction_files)} prediction files to process")
    
    # First, get all sample names
    sample_names = [Path(f).stem for f in prediction_files]
    
    # Process files for each cell type
    for cell_type, col_name in cell_type_cols.items():
        print(f"Processing cell type: {cell_type}")
        
        # Create a dictionary to store all data for this cell type
        cell_type_data = {}
        
        # Process each input file
        for file_path in prediction_files:
            sample_name = Path(file_path).stem
            print(f"Processing {sample_name} for {cell_type}")
            
            try:
                # Read the file with pandas instead of cuDF
                df = pd.read_csv(file_path)
                
                # Get gene names and predictions
                if 'gene_name' in df.columns:
                    genes = df['gene_name'].values
                else:
                    # If there's no gene_name column, check if genes are in the first column
                    first_col = df.columns[0]
                    genes = df[first_col].values
                    # Remove the first column from further processing if it contains gene names
                    df = df.drop(columns=[first_col])
                
                # Check if the cell type column exists
                if col_name in df.columns:
                    predictions = df[col_name].values
                    
                    # Store in dictionary
                    for gene, pred in zip(genes, predictions):
                        if gene not in cell_type_data:
                            cell_type_data[gene] = {}
                        cell_type_data[gene][sample_name] = pred
                else:
                    print(f"Warning: Column '{col_name}' not found in {file_path}")
                
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                continue
        
        # Convert the nested dictionary to a pandas DataFrame
        if cell_type_data:
            result_df = pd.DataFrame.from_dict(cell_type_data, orient='index')
            
            # Ensure all samples are present (fill missing with NaN)
            for sample in sample_names:
                if sample not in result_df.columns:
                    result_df[sample] = np.nan
            
            # Sort columns (samples) alphabetically
            result_df = result_df.sort_index(axis=1)
            
            # Add gene_name column by resetting the index
            result_df.reset_index(inplace=True)
            result_df.rename(columns={'index': 'gene_name'}, inplace=True)
            
            # Save to file
            output_file = os.path.join(output_dir, f"{cell_type}_predictions.csv")
            result_df.to_csv(output_file, index=False)
            print(f"Saved {cell_type} predictions to {output_file}")
        else:
            print(f"No data found for cell type: {cell_type}")

def main():
    # Define directories 
    input_dir = "/fs/ess/PAS2598/h5/3_Predicted_Expressions"
    output_dir = "/fs/ess/PAS2598/h5/4_CellType_Preds"
    
    # Process all predictions
    process_predictions(input_dir, output_dir)
    print("Done!")

if __name__ == "__main__":
    main()