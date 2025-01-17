import h5py
import torch
import numpy as np
import pandas as pd
from pathlib import Path
import re
from scipy.stats import rankdata

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

def load_region_gene_mapping(mapping_path):
    """Load the mapping between H5 regions and genes"""
    print("Loading region-gene mapping...")
    # Ensure gene_id is read as string
    mapping_df = pd.read_csv(mapping_path, dtype={'gene_id': str})

    # Prioritize 'within' relationships
    priority_mapping = mapping_df[mapping_df['relationship'] == 'within']

    # For genes without 'within' relationships, use the closest up/downstream
    remaining_genes = set(mapping_df['gene_id'].astype(str)) - set(priority_mapping['gene_id'].astype(str))
    if remaining_genes:
        remaining_df = mapping_df[mapping_df['gene_id'].isin(remaining_genes)]
        closest_positions = remaining_df.sort_values('distance').groupby('gene_id').first()
        priority_mapping = pd.concat([priority_mapping, closest_positions])

    # Create dictionary mapping H5 keys to gene IDs, ensuring string type
    region_to_gene = dict(zip(priority_mapping['h5_key'], priority_mapping['gene_id'].astype(str)))
    gene_to_region = dict(zip(priority_mapping['gene_id'].astype(str), priority_mapping['h5_key']))

    print(f"Found mappings for {len(gene_to_region)} genes")
    return region_to_gene, gene_to_region

def normalize_predictions(predictions):
    ranks = pd.Series(predictions).rank(method='average', pct=True)
    return ranks.values

def load_enformer_features(h5_path, gene_to_region):
    """Load and process Enformer features from H5 file"""
    features_dict = {}

    with h5py.File(h5_path, 'r') as f:
        print("\nReading features from H5 file...")
        print(f"Found {len(f.keys())} regions in H5 file")

        # Process each gene
        for gene_id, region_key in gene_to_region.items():
            if region_key in f:
                # Get the 4x5313 matrix and average it
                region_features = f[region_key][:]
                avg_features = region_features.mean(axis=0)
                features_dict[str(gene_id)] = avg_features

    print(f"\nLoaded features for {len(features_dict)} genes")

    # Convert to array format, ensuring all keys are strings
    gene_ids = sorted(str(key) for key in features_dict.keys())
    features_array = np.stack([features_dict[gene_id] for gene_id in gene_ids])

    return features_array, gene_ids

def predict_expression(model_path, features):
    """Generate predictions using a trained ctPred model"""
    model = CtPred()
    model.load_state_dict(torch.load(model_path, map_location='cpu'))
    model.eval()

    features_tensor = torch.FloatTensor(features)

    with torch.no_grad():
        predictions = model(features_tensor)

    # Convert predictions to numpy and flatten
    predictions = predictions.numpy().flatten()

    # Normalize predictions to ranks between 0 and 1
    normalized_predictions = normalize_predictions(predictions)

    return normalized_predictions

def main(h5_path, mapping_path, models_dir, output_dir):
    """Main function to process data and generate predictions"""
    # Convert paths to absolute paths
    h5_path = Path(h5_path).resolve()
    mapping_path = Path(mapping_path).resolve()
    models_dir = Path(models_dir).resolve()
    output_dir = Path(output_dir).resolve()

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load region-gene mapping
    region_to_gene, gene_to_region = load_region_gene_mapping(mapping_path)

    # Load features
    print("Loading Enformer features...")
    features, gene_ids = load_enformer_features(h5_path, gene_to_region)

    print(f"\nProcessing predictions for {len(gene_ids)} genes")

    # Process each model
    for model_path in Path(models_dir).glob('*.pt'):
        cell_type = model_path.stem
        print(f"\nProcessing {cell_type}...")

        # Generate predictions
        predictions = predict_expression(str(model_path), features)

        # Create output dataframe
        results_df = pd.DataFrame({
            'gene_id': gene_ids,
            'predicted_expression': predictions
        })

        # Save predictions
        output_path = Path(output_dir) / f"{cell_type}_predictions.csv"
        results_df.to_csv(output_path, index=False)
        print(f"Saved predictions to {output_path}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Generate cell-type specific expression predictions')
    parser.add_argument('--h5_path', required=True, help='Path to Enformer H5 file')
    parser.add_argument('--mapping_path', required=True, help='Path to region-gene mapping CSV')
    parser.add_argument('--models_dir', required=True, help='Directory containing trained models')
    parser.add_argument('--output_dir', required=True, help='Directory to save predictions')

    args = parser.parse_args()
    main(args.h5_path, args.mapping_path, args.models_dir, args.output_dir)
