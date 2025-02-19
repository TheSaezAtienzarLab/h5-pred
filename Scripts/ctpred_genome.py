import torch
import pandas as pd
import sys
import os
from pathlib import Path
import glob
import cupy as cp
sys.path.append('/users/PAS2598/duarte63/GitHub/h5-pred/ctPred')
from ctPred_utils import ctPred

print("Script started")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"CUDA device count: {torch.cuda.device_count()}")
    print(f"Current CUDA device: {torch.cuda.current_device()}")

def load_model(filepath):
    print(f"Loading model from {filepath}")
    model = ctPred()
    # Load the state dict with a map location to handle device mismatch
    state_dict = torch.load(filepath, map_location='cuda')
    model.load_state_dict(state_dict)
    model = model.cuda()  # Move model to GPU
    model.eval()
    return model

def process_model(model_path, features_df):
    try:
        print(f"Processing {model_path.name}")
        model = load_model(model_path)
        model.eval()  # Ensure model is in evaluation mode
        
        # Extract numeric features
        numeric_features = features_df.iloc[:, 1:5314].values
        print(f"Features shape before normalization: {numeric_features.shape}")
        
        # Convert to tensor and move to GPU
        features = torch.tensor(numeric_features, dtype=torch.float32).cuda()
        
        # Normalize features using the same approach as in training
        features_mean = features.mean(dim=0)
        features_std = features.std(dim=0)
        normalized_features = (features - features_mean) / features_std
        print(f"Normalized features shape: {normalized_features.shape}")
        
        # Get predictions
        with torch.no_grad():
            predictions = model(normalized_features)
            print(f"Raw predictions shape: {predictions.shape}")
            
            # Ensure predictions are the right shape (should be [n_samples, 1])
            if len(predictions.shape) != 2 or predictions.shape[1] != 1:
                predictions = predictions.view(-1, 1)
                print(f"Reshaped predictions: {predictions.shape}")
        
        # Convert to CuPy array
        predictions_cp = cp.asarray(predictions.contiguous())
        
        # Create results DataFrame
        cell_type = model_path.stem.replace('_ctPred', '')
        results_df = pd.DataFrame({
            'gene_name': features_df['gene_name'],
            f'{cell_type}_mean_expression_predicted_expression': cp.asnumpy(predictions_cp).flatten()
        })
        
        # Verify output
        n_input_genes = len(features_df['gene_name'].unique())
        n_output_predictions = len(results_df)
        if n_input_genes != n_output_predictions:
            raise ValueError(f"Mismatch in predictions: {n_input_genes} input genes vs {n_output_predictions} predictions")
            
        return results_df
        
    except Exception as e:
        print(f"Error in process_model: {str(e)}")
        raise

def process_sample(features_path, model_dir, output_dir):
    try:
        sample_name = Path(features_path).stem.replace('_processed_mean_features', '')
        print(f"\nProcessing sample: {sample_name}")
        
        print("Loading features...")
        features_df = pd.read_csv(features_path, compression='gzip')
        print(f"Input DataFrame shape: {features_df.shape}")
        
        model_files = list(Path(model_dir).glob('*_ctPred.pt'))
        print(f"Found {len(model_files)} models")
        
        if not model_files:
            print(f"No models found in {model_dir}")
            return False
            
        results = process_model(model_files[0], features_df)
        
        for model_path in model_files[1:]:
            model_results = process_model(model_path, features_df)
            results = pd.merge(results, model_results, on='gene_name')
        
        output_path = os.path.join(output_dir, f"{sample_name}_predicted_expressions.csv.gz")
        print(f"Saving results to {output_path}")
        results.to_csv(output_path, index=False, compression='gzip')
        return True
    
    except Exception as e:
        print(f"Error processing {features_path}: {str(e)}")
        return False

def main():
    print("Starting main function")
    features_dir = "/fs/ess/PAS2598/h5/2_Gene_Features"
    model_dir = "/users/PAS2598/duarte63/GitHub/h5-pred/Model"
    output_dir = "/fs/ess/PAS2598/h5/3_predicted_Expressions"
    
    print(f"Checking directories:")
    print(f"Features dir exists: {os.path.exists(features_dir)}")
    print(f"Model dir exists: {os.path.exists(model_dir)}")
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    feature_files = glob.glob(os.path.join(features_dir, "*_processed_mean_features.csv.gz"))
    print(f"Found {len(feature_files)} samples to process")
    
    if not feature_files:
        print(f"No files found matching pattern in {features_dir}")
        print("Directory contents:")
        print(os.listdir(features_dir))
        return
        
    successful = 0
    failed = 0
    
    for feature_file in feature_files:
        if process_sample(feature_file, model_dir, output_dir):
            successful += 1
        else:
            failed += 1
    
    print("\nProcessing Summary:")
    print(f"Successfully processed: {successful} samples")
    print(f"Failed to process: {failed} samples")

if __name__ == "__main__":
    try:
        print("Starting script execution")
        main()
    except Exception as e:
        print(f"Fatal error: {str(e)}")