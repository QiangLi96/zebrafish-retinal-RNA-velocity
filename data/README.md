# Data

This directory provides information on the datasets used in this study.

Due to their file sizes, the data files are hosted on Google Drive and can be accessed from:

[Download data from Google Drive](https://drive.google.com/drive/folders/1pco-DO5XUinvAA4MROjInsY_Yga90cAP)

## Data files

### `H44_H48_H60_H72_DL_1208_pres.h5ad`

Single-cell transcriptomic dataset used as the input for data preprocessing and subsequent RNA-velocity analysis.

### `H44DL1208_pres_processed_data.h5ad`

Processed dataset generated from `H44_H48_H60_H72_DL_1208_pres.h5ad` using:

[`H44DL1208pres_data_process_1.py`](../H44DL1208pres_data_process_1.py)

This processed dataset is used for downstream vector-field reconstruction and dynamical analyses.
