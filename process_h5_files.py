import h5py
import torch
import numpy as np
from pathlib import Path
from collections import OrderedDict
from tqdm import tqdm
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

class CtPredModel(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.net = torch.nn.ModuleDict({
            '0': torch.nn.Linear(5313, 64),
            '3': torch.nn.Linear(64, 64),
            '6': torch.nn.Linear(64, 64),
            '9': torch.nn.Linear(64, 64),
            '12': torch.nn.Linear(64, 1)
        })

    def forward(self, x):
        x = torch.relu(self.net['0'](x))
        x = torch.relu(self.net['3'](x))
        x = torch.relu(self.net['6'](x))
        x = torch.relu(self.net['9'](x))
        x = self.net['12'](x)
        return x

def load_model(model_path, device='cuda'):
    if not torch.cuda.is_available() and device == 'cuda':
        print("CUDA is not available, falling back to CPU")
        device = 'cpu'
    
    # Load the state dict with map_location to handle GPU/CPU transitions properly
    state_dict = torch.load(model_path, map_location=device)
    
    # Create model and move to device
    model = CtPredModel().to(device)
    
    # Create new state dict with correct device placement
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        v = v.to(device)
        new_state_dict[k] = v
    
    # Load state dict and set to eval mode
    model.load_state_dict(new_state_dict)
    model.eval()
    
    if device == 'cuda':
        torch.cuda.empty_cache()  # Clear any unused memory
    
    return model

def process_h5_batch(h5_file, model, batch_size=100, device='cuda'):
    results = {}
    with h5py.File(h5_file, 'r') as f:
        datasets = list(f.keys())

        for i in tqdm(range(0, len(datasets), batch_size)):
            batch_datasets = datasets[i:i + batch_size]
            batch_data = []
            batch_names = []

            for dataset_name in batch_datasets:
                data = f[dataset_name][:]
                batch_data.append(data)
                batch_names.append(dataset_name)

            batch_tensor = torch.tensor(np.array(batch_data), dtype=torch.float32).to(device)
            batch_size, num_samples, num_features = batch_tensor.shape
            batch_tensor = batch_tensor.reshape(-1, num_features)

            with torch.no_grad():
                predictions = model(batch_tensor)
                predictions = predictions.reshape(batch_size, num_samples)
                predictions = predictions.cpu().numpy()

            for name, pred in zip(batch_names, predictions):
                results[name] = pred

    return results

def predict_all_samples(h5_dir, models_dir, output_dir, device='cuda'):
    h5_dir = Path(h5_dir)
    output_dir = Path(output_dir)
    h5_files = list(h5_dir.glob('*.h5'))
    model_files = list(Path(models_dir).glob('*_ctPred.pt'))
    
    # Track completed samples
    completed_samples = {d.name for d in output_dir.iterdir() if d.is_dir()} if output_dir.exists() else set()
    logging.info(f"Found {len(h5_files)} H5 files and {len(model_files)} cell type models")
    logging.info(f"Already completed {len(completed_samples)} samples")

    for h5_file in h5_files:
        sample_name = h5_file.stem
        if sample_name in completed_samples:
            logging.info(f"Skipping completed sample: {sample_name}")
            continue
            
        sample_output_dir = output_dir / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        logging.info(f"\nProcessing sample: {sample_name}")

        for model_path in model_files:
            cell_type = model_path.stem.replace('_mean_expression_ctPred', '')
            logging.info(f"Processing cell type: {cell_type}")

            try:
                model = load_model(model_path, device)
                results = process_h5_batch(h5_file=h5_file, model=model, device=device)

                output_file = sample_output_dir / f"{cell_type}_predictions.h5"
                with h5py.File(output_file, 'w') as f:
                    for region_name, predictions in results.items():
                        f.create_dataset(region_name, data=predictions)

                logging.info(f"Saved predictions to {output_file}")

            except Exception as e:
                logging.error(f"Error processing {cell_type} for sample {sample_name}: {str(e)}")
                continue

            finally:
                del model
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

def main():
    h5_dir = "/fs/ess/PAS2598/h5/Enformer_output_4bins"
    models_dir = "/users/PAS2598/duarte63/GitHub/h5-pred/Model"
    output_dir = "predictions"
    
    # Use CUDA if available, else CPU
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    predict_all_samples(h5_dir, models_dir, output_dir, device)

if __name__ == "__main__":
    main()