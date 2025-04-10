# In Silico Vaccine Design Automation

A NextFlow-based pipeline for automated computational vaccine design targeting avian influenza viruses, based on epitope prediction and molecular dynamics simulations.

## Overview

This pipeline automates the process of designing epitope-based vaccines for H5N1 avian influenza virus using computational methods. The workflow implements the methodology described in:

> Tambunan et al. Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions. *Bioinformatics and Biology Insights* 2016:10 27-35.

The pipeline incorporates the following steps:
1. Retrieve protein sequences from NCBI
2. Predict B-cell epitopes
3. Predict T-cell epitopes (MHC class I and II binding)
4. Select high-quality epitopes based on scoring
5. Design a multi-epitope vaccine construct with appropriate linkers
6. Evaluate vaccine properties (optional molecular dynamics simulation)

## Requirements

- Nextflow (>=21.10)
- Container technology:
  - Docker or Singularity
- R (>=4.0) packages:
  - Biostrings
  - rentrez
  - seqinr
  - httr
- Python (>=3.8) packages:
  - biopython
  - pandas
  - numpy
  - requests

## Quick Start

```bash
# Clone this repository
git clone https://github.com/username/in-silico-vaccine-design.git
cd in-silico-vaccine-design

# Run with default parameters
nextflow run main.nf

# Run with custom protein accession
nextflow run main.nf --accession "BAL61222.1"

# Run on a SLURM cluster
nextflow run main.nf -profile slurm
```

## Pipeline Parameters

Edit the `params.yaml` file or provide parameters at runtime:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `accession` | NCBI protein accession number | BAL61222.1 |
| `outdir` | Output directory | results |
| `bcell_method` | B-cell prediction method | Bepipred |
| `tcell_mhci_method` | MHC-I prediction method | NetMHCpan |
| `tcell_mhcii_method` | MHC-II prediction method | NetMHCIIpan |
| `run_md` | Run molecular dynamics simulation | false |

## Output

The pipeline generates structured output in the specified output directory:

- `sequences/` - Retrieved protein sequences in FASTA format
- `epitopes/` - Predicted epitopes for B-cells, T-cells (MHC-I and MHC-II)
- `vaccine_construct.fasta` - Final designed vaccine sequence
- `evaluation/` - Evaluation reports for the vaccine construct
- `molecular_dynamics/` - Molecular dynamics simulation results (if enabled)

## Customization

You can customize HLA alleles for T-cell epitope prediction and add additional prediction methods by editing the configuration file:

```bash
# Example of running with custom HLA alleles
nextflow run main.nf --mhci_alleles "HLA-A*02:01,HLA-B*07:02,HLA-B*35:01"
```

## Citation

If you use this pipeline, please cite:

```
Tambunan et al. Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions. Bioinformatics and Biology Insights 2016:10 27-35.
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.