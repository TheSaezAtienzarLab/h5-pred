import torch
import pandas as pd
import sys
import os
from pathlib import Path
import glob
import cupy as cp
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('prediction_processing.log'),
        logging.StreamHandler()
    ]
)

# Print CUDA information
logging.info(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    logging.info(f"CUDA device count: {torch.cuda.device_count()}")
    logging.info(f"Current CUDA device: {torch.cuda.current_device()}")

def get_cell_type_columns():
    """Return dictionary of cell types and their corresponding column names."""
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

def load_model(filepath):
    """Load a PyTorch model from file."""
    logging.info(f"Loading model from {filepath}")
    
    # Import ctPred class
    sys.path.append('/users/PAS2598/duarte63/GitHub/h5-pred/ctPred')
    from ctPred_utils import ctPred
    
    model = ctPred()
    # Load the state dict with map location to handle device mismatch
    state_dict = torch.load(filepath, map_location='cuda')
    model.load_state_dict(state_dict)
    model = model.cuda()  # Move model to GPU
    model.eval()  # Set model to evaluation mode
    return model

def process_model(model_path, features_df):
    """Process a single model's predictions."""
    try:
        logging.info(f"Processing {model_path.name}")
        model = load_model(model_path)
        
        # Extract numeric features
        numeric_features = features_df.iloc[:, 1:5314].values
        logging.info(f"Features shape before normalization: {numeric_features.shape}")
        
        # Convert to tensor and move to GPU
        features = torch.tensor(numeric_features, dtype=torch.float32).cuda()
        
        # Normalize features
        features_mean = features.mean(dim=0)
        features_std = features.std(dim=0)
        normalized_features = (features - features_mean) / (features_std + 1e-8)
        logging.info(f"Normalized features shape: {normalized_features.shape}")
        
        # Get predictions
        with torch.no_grad():
            predictions = model(normalized_features)
            logging.info(f"Raw predictions shape: {predictions.shape}")
            
            # Ensure predictions are the right shape
            if len(predictions.shape) != 2 or predictions.shape[1] != 1:
                predictions = predictions.view(-1, 1)
                logging.info(f"Reshaped predictions: {predictions.shape}")
        
        # Convert to CuPy array
        predictions_cp = cp.asarray(predictions.contiguous())
        
        # Create results DataFrame
        cell_type = model_path.stem.replace('_ctPred', '')
        results_df = pd.DataFrame({
            'gene_name': features_df['gene_name'],
            f'{cell_type}_mean_expression_predicted_expression': cp.asnumpy(predictions_cp).flatten()
        })
        
        # Verify output
        if len(results_df) != len(features_df):
            raise ValueError(f"Mismatch in predictions: {len(features_df)} input rows vs {len(results_df)} predictions")
            
        return results_df
        
    except Exception as e:
        logging.error(f"Error in process_model: {str(e)}")
        raise

def process_sample(features_path, model_dir, output_dir):
    """Process a single sample through all models."""
    try:
        sample_name = Path(features_path).stem.replace('_processed_mean_features', '')
        logging.info(f"\nProcessing sample: {sample_name}")
        
        # Load and preprocess features
        logging.info("Loading features...")
        features_df = pd.read_csv(features_path, compression='gzip')
        logging.info(f"Input DataFrame shape: {features_df.shape}")
        
        # Check for and handle duplicates
        duplicate_genes = features_df[features_df['gene_name'].duplicated()]['gene_name'].unique()
        if len(duplicate_genes) > 0:
            logging.warning(f"Found {len(duplicate_genes)} duplicate genes in input. Keeping first occurrence.")
            if len(duplicate_genes) < 10:
                logging.info(f"Duplicate genes: {', '.join(duplicate_genes)}")
            features_df = features_df.drop_duplicates(subset='gene_name', keep='first')
            logging.info(f"DataFrame shape after removing duplicates: {features_df.shape}")
        
        # Process with each model
        model_files = list(Path(model_dir).glob('*_ctPred.pt'))
        logging.info(f"Found {len(model_files)} models")
        
        if not model_files:
            logging.error(f"No models found in {model_dir}")
            return False
            
        # Process first model
        results = process_model(model_files[0], features_df)
        
        # Process remaining models
        for model_path in model_files[1:]:
            model_results = process_model(model_path, features_df)
            results = pd.merge(results, model_results, on='gene_name')
        
        # Save results
        output_path = os.path.join(output_dir, f"{sample_name}_predicted_expressions.csv.gz")
        logging.info(f"Saving results to {output_path}")
        results.to_csv(output_path, index=False, compression='gzip')
        
        # Verify final output
        logging.info(f"Final output shape: {results.shape}")
        return True
    
    except Exception as e:
        logging.error(f"Error processing {features_path}: {str(e)}")
        return False

def output_exists(sample_name, output_dir):
    """Check if the output file already exists for a given sample."""
    output_path = os.path.join(output_dir, f"{sample_name}_predicted_expressions.csv.gz")
    return os.path.exists(output_path)

def main():
    """Main function to process all samples."""
    logging.info("Starting main function")
    
    # Set up directories
    features_dir = "/fs/ess/PAS2598/h5/2_Gene_Features"
    model_dir = "/users/PAS2598/duarte63/GitHub/h5-pred/Model"
    output_dir = "/fs/ess/PAS2598/h5/3_Predicted_Expressions"
    
    # Verify directories exist
    logging.info("Checking directories:")
    logging.info(f"Features dir exists: {os.path.exists(features_dir)}")
    logging.info(f"Model dir exists: {os.path.exists(model_dir)}")
    
    # Create output directory if needed
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Get list of files to process
    feature_files = glob.glob(os.path.join(features_dir, "*_processed_mean_features.csv.gz"))
    logging.info(f"Found {len(feature_files)} samples to process")
    
    if not feature_files:
        logging.error(f"No files found matching pattern in {features_dir}")
        logging.info("Directory contents:")
        logging.info(os.listdir(features_dir))
        return
    
    # Process all samples
    successful = 0
    failed = 0
    skipped = 0
    
    for feature_file in feature_files:
        sample_name = Path(feature_file).stem.replace('_processed_mean_features', '')
        
        # Check if output already exists
        if output_exists(sample_name, output_dir):
            logging.info(f"Skipping {sample_name} - output already exists")
            skipped += 1
            continue
            
        if process_sample(feature_file, model_dir, output_dir):
            successful += 1
        else:
            failed += 1
    
    # Print summary
    logging.info("\nProcessing Summary:")
    logging.info(f"Successfully processed: {successful} samples")
    logging.info(f"Failed to process: {failed} samples")
    logging.info(f"Skipped (already existed): {skipped} samples")

if __name__ == "__main__":
    try:
        logging.info("Starting script execution")
        main()
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")