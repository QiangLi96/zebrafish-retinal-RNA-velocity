import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dynamo as dyn
import os
import anndata as Anndata

# type set
from typing import Dict, List, Any, Tuple
from pandas import DataFrame

dyn.dynamo_logger.main_silence()


###
# ----------------------------------------------------------------------------------------------

## Data Cleaning
adata: Anndata = dyn.read("zebrafish_H44atoh7_1/data/H44_H48_H60_H72_DL_1208_pres.h5ad")
# umap1: DataFrame = adata.obs['umap_1']
# umap2: DataFrame = adata.obs['umap_2']
adata.layers['spliced'] = adata.layers['spliced'].astype('int64')
adata.layers['unspliced'] = adata.layers['unspliced'].astype('int64')
adata.X = adata.layers['spliced']
adata.X = adata.X.astype('float32')
# adata.obs['timebatch'] = adata.obs['timebatch'].astype('str')
# adata.obs['timebatch'] = adata.obs['timebatch'].astype('category')
del adata.layers['ambiguous']
del adata.layers['matrix']


## dynamo recipe data
dyn.pp.recipe_monocle(adata)
dyn.tl.dynamics(adata, model='stochastic', cores=3)
# or dyn.tl.dynamics(adata, model='deterministic')
# or dyn.tl.dynamics(adata, model='stochastic', est_method='negbin')

## dynamo reconstruct vector feild in umap space
dyn.tl.reduceDimension(adata)
adata.obsm['X_umap'][:, 0] = adata.obs['umap_1']
adata.obsm['X_umap'][:, 1] = adata.obs['umap_2']

# adata.obsm['X_umap'] = np.load('zebrafish_H44atoh7_1\Figures_H44_DL_1208_pres\H44_DL_1208_pres_umap_auto_1.npy')

dyn.pl.umap(adata, color='Cell_type', show_legend='on data')


#dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'})
# dyn.tl.cell_velocities(adata, method='cosine', other_kernels_dict={'transform': 'sqrt'})
# dyn.tl.cell_velocities(adata, method='fp', other_kernels_dict={'transform': 'sqrt'})

# You can check the confidence of cell-wise velocity to understand how reliable the recovered velocity is across cells or even correct velocty based on some prior:
dyn.tl.cell_wise_confidence(adata)
dyn.tl.confident_cell_velocities(adata, group='Cell_type', lineage_dict={'pre1': ['RGCs', 'GABAergic ACs', 'PRs'], 'PRpre': 'PRs', 'pre2': ['Glycinergic ACs', 'BCs'], 'HCpre': 'HCs'})

## dynamo reconstruct vector feild in umap space
dyn.pl.cell_wise_vectors(adata, color=['Cell_type'], basis='umap', show_legend='on data', quiver_length=6,quiver_size=6, pointsize=0.1, show_arrowed_spines=False)

dyn.pl.streamline_plot(adata, color=['Cell_type'], basis='umap', show_legend='on data',
            show_arrowed_spines=True)
# you can set `verbose = 1/2/3` to obtain different levels of running information of vector field reconstruction
# M is the number of control points
dyn.vf.VectorField(adata, basis='umap', n=100, M=1000,pot_curl_div=True)
dyn.pl.topography(adata, basis='umap', background='white', color=['ntr', 'Cell_type'],streamline_color='black',
                  show_legend='on data', frontier=True,save_show_or_return='return')
# black: absorbing fixed points;
# red: emitting fixed points;
# blue: unstable fixed points.

adata = dyn.read_h5ad("zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")
## dynamo reconstruct vector feild in pca space
dyn.tl.cell_velocities(adata, basis='pca', method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.pl.cell_wise_vectors(adata, color=['Cell_type'], basis='pca', show_legend='on data', quiver_length=6, quiver_size=6, pointsize=0.1, show_arrowed_spines=False)
dyn.pl.streamline_plot(adata, color=['Cell_type'], basis='pca', show_legend='on data',
            show_arrowed_spines=True)
# you can set `verbose = 1/2/3` to obtain different levels of running information of vector field reconstruction
# M is the number of control points
dyn.vf.VectorField(adata, basis='pca',n=100,M=1000, pot_curl_div=True,method='SparseVFC')
dyn.pl.topography(adata, basis='pca',background='white', color=['ntr', 'Cell_type'],streamline_color='black', show_legend='on data', frontier=True)


# ----------------------------------------------------------------------------------------------
###

### Dynamo save utility
# there may be intermediate results stored in adata.uns that can may lead to errors when writing the h5ad object.
# call dyn.cleanup(adata) first to remove these data objects before saving the adata object.
dyn.cleanup(adata)
# # call AnnData write_h5ad to save the entire adata information.
adata.write_h5ad("zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")
# ----------------------------------------------------------------------------------------------