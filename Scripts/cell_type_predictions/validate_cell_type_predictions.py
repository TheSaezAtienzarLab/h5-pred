import pandas as pd
from pathlib import Path
import glob
import os

def validate_output_files(output_dir):
    """Validate the structure of output files and print summary statistics."""
    print("\nValidating output files...")
    
    for file_path in glob.glob(os.path.join(output_dir, "*_predictions.csv")):
        cell_type = Path(file_path).stem.replace('_predictions', '')
        
        # Read the first few rows of the file
        df = pd.read_csv(file_path, nrows=5)
        
        # Get structure information
        n_samples = len(df.columns) - 1  # -1 for gene_name column
        n_genes = sum(1 for _ in open(file_path)) - 1  # -1 for header
        
        print(f"\nFile: {Path(file_path).name}")
        print(f"Structure:")
        print(f"- Number of samples (columns): {n_samples}")
        print(f"- Number of genes (rows): {n_genes}")
        print("\nFirst few rows:")
        print(df.head())
        print("\nColumn names:")
        print(df.columns.tolist())
        
        # Basic data validation
        print("\nValidation checks:")
        print(f"- Has 'gene_name' column: {'gene_name' in df.columns}")
        print(f"- All numeric values (excluding gene_name): {df.drop('gene_name', axis=1).apply(pd.to_numeric, errors='coerce').notnull().all().all()}")
        print(f"- No duplicate columns: {len(df.columns) == len(set(df.columns))}")
        print(f"- No missing gene names: {not df['gene_name'].isna().any()}")
        print(f"- No duplicate gene names: {not df['gene_name'].duplicated().any()}")
        print("-" * 80)

def main():
    output_dir = "/fs/ess/PAS2598/h5/4_CellType_Preds"
    validate_output_files(output_dir)

if __name__ == "__main__":
    main() 