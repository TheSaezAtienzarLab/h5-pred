import h5py
import torch
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

class CtPred(torch.nn.Module):
    def __init__(self, input_size=5313):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_size, 64),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.05),
            torch.nn.Linear(64, 64),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.05),
            torch.nn.Linear(64, 64),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.05),
            torch.nn.Linear(64, 64),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.05),
            torch.nn.Linear(64, 1)
        )

    def forward(self, x):
        return self.net(x)

def load_gene_annotations(anno_path):
    """Load gene annotations from tab-delimited file"""
    print("Loading gene annotations...")
    anno_df = pd.read_csv(anno_path, sep='\t')
    # Create a dictionary mapping gene_ids to their genomic coordinates
    gene_coords = anno_df.set_index('gene_id').to_dict('index')
    print(f"Loaded {len(gene_coords)} gene annotations")
    return gene_coords

def get_region_key(chrom, start, end):
    """Generate H5 key format for a genomic region"""
    return f"chr{chrom}:{start}-{end}"

def load_enformer_features(h5_path, gene_coords):
    """Load and process Enformer features from a single H5 file"""
    features_dict = {}
    
    try:
        with h5py.File(h5_path, 'r') as f:
            for gene_id, coords in gene_coords.items():
                region_key = get_region_key(coords['chr'], coords['start'], coords['end'])
                if region_key in f:
                    region_features = f[region_key][:]
                    avg_features = region_features.mean(axis=0)  # Average the 4 bins
                    features_dict[gene_id] = avg_features
    except Exception as e:
        print(f"Error processing {h5_path}: {str(e)}")
        return None, None

    # Convert to array format
    gene_ids = sorted(features_dict.keys())
    if gene_ids:  # Only if we found any features
        features_array = np.stack([features_dict[gene_id] for gene_id in gene_ids])
        return features_array, gene_ids
    return None, None

def normalize_predictions(predictions):
    """Normalize predictions to [0,1] range"""
    return pd.Series(predictions).rank(method='average', pct=True).values

def predict_expression(model_path, features):
    """Generate predictions using a trained ctPred model"""
    try:
        model = CtPred()
        model.load_state_dict(torch.load(model_path, map_location='cpu'))
        model.eval()

        features_tensor = torch.FloatTensor(features)
        
        with torch.no_grad():
            predictions = model(features_tensor)
        
        predictions = predictions.numpy().flatten()
        normalized_predictions = normalize_predictions(predictions)
        
        return normalized_predictions
    except Exception as e:
        print(f"Error in prediction for model {model_path}: {str(e)}")
        return None

def process_sample(h5_path, gene_coords, models_dir, output_dir):
    """Process a single sample for all cell types"""
    sample_id = h5_path.stem
    
    # Create sample-specific output directory
    sample_output_dir = Path(output_dir) / sample_id
    sample_output_dir.mkdir(parents=True, exist_ok=True)

    # Load features for this sample
    features, gene_ids = load_enformer_features(h5_path, gene_coords)
    
    if features is None:
        print(f"No features found for {sample_id}, skipping...")
        return

    # Process each model for this sample
    for model_path in Path(models_dir).glob('*_ctPred.pt'):
        cell_type = model_path.stem.replace('_mean_expression_ctPred', '')
        
        # Generate predictions
        predictions = predict_expression(str(model_path), features)
        
        if predictions is not None:
            # Create output dataframe
            results_df = pd.DataFrame({
                'gene_id': gene_ids,
                'predicted_expression': predictions
            })

            # Save predictions
            output_path = sample_output_dir / f"{cell_type}_predictions.csv"
            results_df.to_csv(output_path, index=False)
            print(f"Saved predictions for {sample_id} - {cell_type}")

def main():
    # Define paths
    h5_dir = "/fs/ess/PAS2598/h5/Enformer_output_4bins"
    models_dir = "/users/PAS2598/duarte63/GitHub/h5-pred/Model"
    anno_path = "/users/PAS2598/duarte63/GitHub/h5-pred/gene_anno.txt"
    output_dir = "predictions"

    # Create main output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load gene annotations
    gene_coords = load_gene_annotations(anno_path)

    # Get list of all H5 files
    h5_files = list(Path(h5_dir).glob('*.h5'))
    print(f"\nFound {len(h5_files)} H5 files")

    # Process each H5 file
    for h5_path in tqdm(h5_files, desc="Processing samples"):
        process_sample(h5_path, gene_coords, models_dir, output_dir)

if __name__ == "__main__":
    main()