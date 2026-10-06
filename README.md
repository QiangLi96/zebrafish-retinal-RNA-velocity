# Dynamical landscape reconstruction of retinal differentiation from stabilized RNA-velocity vector fields

This repository contains the Python code and data-processing workflow associated with the study:

**Dynamical landscape reconstruction of retinal differentiation from stabilized RNA-velocity vector fields**

The computational framework reconstructs a continuous RNA-velocity vector field of zebrafish retinal development and performs downstream dynamical analyses, including fixed-point identification, Jacobian-based regulatory analysis, least action path (LAP) analysis, gene-velocity synchronization, and in silico gene perturbation.

## Data

The datasets used in this study are provided through Google Drive because of their file sizes:

[Download data from Google Drive](https://drive.google.com/drive/folders/1pco-DO5XUinvAA4MROjInsY_Yga90cAP)

The `data/` directory contains additional information on the datasets and preprocessing procedure.

The main data files are:

- `H44_H48_H60_H72_DL_1208_pres.h5ad`
- `H44DL1208_pres_processed_data.h5ad`

The latter is generated from `H44_H48_H60_H72_DL_1208_pres.h5ad` using:

`H44DL1208pres_data_process_1.py`

## Software dependency

The analyses are based in part on the Python package [dynamo](https://github.com/aristoteleo/dynamo-release).

Before running the vector-field analyses, install `dynamo-release` following its official instructions.

### Important: modified `scVectorField.py`

To reproduce the stabilized vector-field reconstruction used in this study, the original

```text
dynamo/vectorfield/scVectorField.py

in the installed `dynamo-release` package should be replaced with the modified version provided in this repository:
`scVectorField.py`
