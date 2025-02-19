# h5-pred: GPU-Accelerated Genomic Feature Processing

## Overview
A collection of Python scripts for efficient processing and reorganization of genomic feature data using GPU acceleration. This repository focuses on handling gene region features, H5 files, and cell type predictions with CUDA-enabled processing via cupy and cuDF.

## Scripts

### 1. `calculate_gene_region_means.py`
Calculates mean features for gene regions using GPU acceleration.
- Processes gene windows and their associated features
- Groups regions by gene
- Outputs averaged features per gene
- Uses: cupy, pandas, numpy, tqdm

### 2. `gpu_h5_to_csv_means.py`
Converts H5 files containing gene matrices to CSV format with calculated means.
- Processes H5 files with 4x5313 matrices per gene
- Calculates mean features using GPU
- Exports results to CSV with gene metadata
- Uses: h5py, cupy, pandas

### 3. `split_predictions_by_cell_type.py`
Reorganizes prediction data by cell type using GPU acceleration.
- Splits combined prediction files into cell-type-specific files
- Efficiently processes large datasets using GPU
- Uses: cuDF, cupy

## Requirements
- CUDA-capable GPU
- Python 3.x
- Dependencies:
  - cupy
  - cuDF
  - pandas
  - numpy
  - h5py
  - tqdm

## Installation

1. Clone the repository:

```
git clone https://github.com/username/h5-pred.git
cd h5-pred
```

2. Install required packages:

```
conda install cupy-cuda11x pandas numpy h5py tqdm
```

Note: Replace `cuda11x` with your CUDA version

## Usage
Each script can be run independently from the command line. Some scripts use hardcoded paths that need to be modified before use.

## Directory Structure

h5-pred/
├── calculate_gene_region_means.py
├── gpu_h5_to_csv_means.py
├── split_predictions_by_cell_type.py
└── .gitignore

## Ignored Directories
- `/ctPred`: Directory containing cell type prediction data (ignored in git)

## Contributing
[Add contribution guidelines if applicable]

## License
 