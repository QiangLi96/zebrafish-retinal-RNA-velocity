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

`H44DL1208_pres_processed_data.h5ad` is generated from `H44_H48_H60_H72_DL_1208_pres.h5ad` using:

`H44DL1208pres_data_process_1.py`

The processed dataset is used for subsequent RNA-velocity vector-field reconstruction and downstream dynamical analyses.

## Software dependency

The analyses are based in part on the Python package [dynamo](https://github.com/aristoteleo/dynamo-release).

Before running the vector-field analyses, install `dynamo-release` following its official installation instructions.

### Important: modified `scVectorField.py`

To reproduce the stabilized vector-field reconstruction used in this study, the original

```text
dynamo/vectorfield/scVectorField.py
```

in the installed `dynamo-release` package should be replaced with the modified version provided in this repository:

```text
scVectorField.py
```

The modified implementation introduces the numerical treatment used in this study to stabilize the ill-conditioned kernel coefficient system during RNA-velocity vector-field reconstruction.

After replacement, the analysis scripts can be executed using the modified `dynamo` implementation.

## Analysis workflow

The main computational workflow consists of:

1. Preprocessing of the integrated zebrafish retinal single-cell dataset.
2. RNA-velocity estimation and visualization.
3. SVD-regularized RNA-velocity vector-field reconstruction.
4. Fixed-point identification and dynamical-landscape analysis.
5. Jacobian-based inference of local gene regulatory interactions.
6. Least action path analysis.
7. Analysis of gene-velocity synchronization along LAPs.
8. In silico gene perturbation and trajectory analysis.

## Reproducibility

To reproduce the analyses:

1. Download the datasets from the Google Drive link above.
2. Place the data files in an appropriate local directory.
3. Update the corresponding data paths in the Python scripts if necessary.
4. Install `dynamo-release` and the required Python dependencies.
5. Replace the original

   ```text
   dynamo/vectorfield/scVectorField.py
   ```

   with the `scVectorField.py` provided in this repository.
6. Run `H44DL1208pres_data_process_1.py` if regeneration of the processed dataset is required.
7. Run the corresponding analysis scripts for vector-field reconstruction, fixed-point analysis, Jacobian analysis, least action paths, gene-velocity synchronization, and perturbation analysis.

## Methodological overview

The vector-field reconstruction is based on RNA-velocity estimates derived from single-cell transcriptomic data. A continuous and differentiable vector field is reconstructed in gene-expression state space using an RKHS-based formulation.

To address numerical instability caused by the ill-conditioned kernel coefficient system, singular value decomposition (SVD)-based regularization is introduced during coefficient estimation. The stabilized vector field is subsequently used to characterize developmental trajectories, fixed points, attractor states, local Jacobians, least action paths, and perturbation responses.


## Citation

If you use the code, data-processing workflow, or modified vector-field implementation provided in this repository, please cite the associated study:

> Qiang Li, Mei Wang, Yan Li, Jie He, and Peng Ji.  
> *Dynamical landscape reconstruction of retinal differentiation from stabilized RNA-velocity vector fields.*


