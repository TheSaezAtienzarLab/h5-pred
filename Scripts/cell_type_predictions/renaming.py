import pandas as pd
import os
from pathlib import Path
import glob
from tqdm import tqdm

def get_cell_type_columns():
    """Return dictionary mapping cell types to their full and new column names."""
    cell_types = [
        'pericytes', 'schwann', 'lymphocytes', 'microglia', 'meninges',
        'astrocytes', 'ependymal', 'opc', 'oligodendrocytes', 'neurons', 'endothelial'
    ]
    
    return {
        f"{ct}_mean_expression_mean_expression_predicted_expression": f"{ct}_predicted"
        for ct in cell_types
    }

def process_files(input_dir):
    """Process all prediction files."""
    prediction_files = glob.glob(os.path.join(input_dir, "*_predicted_expressions.csv.gz"))
    print(f"Found {len(prediction_files)} files to process")
    
    cell_type_cols = get_cell_type_columns()
    
    for file_path in tqdm(prediction_files, desc="Processing files"):
        sample_name = Path(file_path).stem.split('.')[0]  # Get name before first .csv
        output_file = os.path.join(os.path.dirname(file_path), f"{sample_name}.csv")
        
        try:
            # Read data
            df = pd.read_csv(file_path, compression='gzip')
            
            # Rename columns
            for old_col, new_col in cell_type_cols.items():
                if old_col in df.columns:
                    df = df.rename(columns={old_col: new_col})
            
            # Write results
            df.to_csv(output_file, index=False)
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue

def main():
    input_dir = "/fs/ess/PAS2598/h5/3_Predicted_Expressions"
    process_files(input_dir)
    print("Done!")

if __name__ == "__main__":
    main()
