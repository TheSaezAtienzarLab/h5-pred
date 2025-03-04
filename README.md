# Description

## scPrediXcan-DataFlow: A Comprehensive Pipeline for Cell-Type Specific TWAS Data Preparation

A modular, GPU-accelerated pipeline for processing multi-omic data through the entire scPrediXcan workflow - from raw genotypes to cell-type specific gene expression predictions. This repository provides the infrastructure to prepare all necessary input files for transcriptome-wide association studies (TWAS) with single-cell resolution, enabling researchers to identify cell-type specific disease mechanisms.
The pipeline integrates genotype processing, Enformer deep learning outputs, and ctPred gene expression predictions into standardized formats ready for association testing. Designed for both high-performance computing environments and standalone workstations, scPrediXcan-DataFlow bridges the gap between raw genomic data and meaningful insights into the cellular basis of complex diseases.

### Key Features:

Complete preparation of all scPrediXcan input files
GPU-optimized processing of large genomic datasets
Cell-type specific gene expression prediction
SLURM integration for high-performance computing
Modular design for flexible workflow adaptation

## Developed at Ohio State University for researchers exploring the genetic basis of disease at single-cell resolution.

# scPrediXcan Pipeline Documentation

This documentation explains the pipeline for generating input files for scPrediXcan, a framework for cell-type-specific transcriptome-wide association studies (TWAS) that leverages single-cell data and deep learning.

## Overview

The pipeline is organized into four major components that process genomic data into the formats required for scPrediXcan:

1. **Genotype Processing**: Prepares genotype data from VCF files
2. **Enformer Processing**: Processes Enformer deep learning model outputs
3. **Cell Type Prediction**: Runs and processes cell-type specific gene expression predictions
4. **scPrediXcan Inputs**: Organizes final files for scPrediXcan linearization step

## Pipeline Components

### 1. Genotype Processing

This component processes raw genotype data from VCF files into the formats required for scPrediXcan.

**Key Scripts:**
- `extract_and_process_genotype.sh`: Extracts genotypes from VCF files, converts to dosage values (0-2)
- `predictdb_genotype_data_final.sh`: Reformats variant IDs (removing "chr" prefix, adding "_GRCh38" suffix)
- `calculate_maf.py/calculate_maf.sh`: Calculates Minor Allele Frequency (MAF) for SNPs
- `rsid_retrieval.py/rsid_retrieval.sh`: Retrieves reference SNP IDs (rsIDs) from dbSNP

**Data Flow:**
```
VCF files → Extract genotypes → Format variant IDs → Calculate MAF → Add rsIDs
```

**Outputs:**
- Genotype files in format required by PredictDb
- Variant annotation files with MAF and rsIDs

### 2. Enformer Processing

This component processes outputs from Enformer, a deep learning model that predicts epigenomic features from DNA sequences.

**Key Scripts:**
- `gpu_h5_to_csv_means.py/gpu_h5_to_csv_means_slurm.sh`: Converts HDF5 files to CSV format
- `calculate_gene_region_means.py/calculate_gene_region_means_slurm.sh`: Calculates mean epigenomic features for gene regions

**Data Flow:**
```
Enformer H5 files → GPU-accelerated processing → Gene region means
```

**Outputs:**
- Processed CSV files with gene-level epigenomic features
- Gene region mean features used as input for cell-type prediction

### 3. Cell Type Prediction

This component runs the ctPred model to predict cell-type-specific gene expression and processes the outputs.

**Key Scripts:**
- `ctpred_genome.py/ctpred_genome_slurm.sh`: Predicts gene expression using genomic features
- `split_predictions_by_cell_type.py/split_predictions_by_cell_type_slurm.sh`: Organizes predictions by cell type
- `renaming.py`: Renames columns for consistency
- `transform_expression_format.py`: Transforms files to PredictDb format
- `validate_cell_type_predictions.py`: Validates output structure and content

**Data Flow:**
```
Gene features → ctPred model → Split by cell type → Rename columns → Format for PredictDb
```

**Outputs:**
- Cell-type-specific gene expression predictions
- Files formatted for use in scPrediXcan linearization step

### 4. scPrediXcan Inputs

This component organizes the final processed files needed for step 2 of scPrediXcan (linearization of ctPred into ℓ-ctPred).

**Directories:**
- `gene_annotation/`: Gene annotation files (Gene_anno.txt)
- `snp_annotation/`: SNP annotation files from genotype processing
- `genotype/`: Final processed genotype files
- `gene_expression/`: Cell-type-specific gene expression predictions

These files serve as inputs to the PredictDb-nextflow pipeline, which creates the SNP-based elastic net models (ℓ-ctPred) used in the final scPrediXcan association tests.

## Connection to scPrediXcan Framework

This pipeline implements the data preparation for the first two steps of the three-step scPrediXcan framework:

1. **Step 1**: Training ctPred (cell-type-specific gene expression prediction model)
   - Uses Enformer as a feature extractor
   - Predicts cell-type-specific gene expression from DNA sequences

2. **Step 2**: Linearizing ctPred into ℓ-ctPred
   - Creates SNP-based elastic net models from ctPred predictions
   - Requires the four types of input files organized in our pipeline

3. **Step 3**: Performing association tests with S-PrediXcan
   - Uses ℓ-ctPred and GWAS summary statistics
   - Identifies genes associated with traits/diseases

## Running the Pipeline

For detailed instructions on how to run each component of the pipeline, refer to the slurm scripts in the respective directories:

1. For genotype processing: `genotype_processing/scripts/`
2. For Enformer processing: `enformer_processing/scripts/`
3. For cell type prediction: `cell_type_predictions/`
4. For running PredictDb: See `README.md` and consult the PredictDb-nextflow documentation

## Technical Considerations

- **GPU Requirements**: Several components (Enformer processing, ctPred) require GPU acceleration
- **Memory Usage**: Processing genomes and large H5 files requires significant memory (32GB+)
- **Storage**: The pipeline generates large intermediate files, especially during Enformer processing
- **SLURM Integration**: All major scripts have SLURM versions for high-performance computing environments

This pipeline enables researchers to leverage the power of deep learning and single-cell data to conduct cell-type-specific transcriptome-wide association studies, providing insights into the cellular mechanisms of complex diseases.
