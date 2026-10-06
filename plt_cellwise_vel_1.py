import warnings
warnings.filterwarnings('ignore')
import dynamo as dyn
import numpy as np
import matplotlib.pyplot as plt

data_path = 'E:/python_work_PyCharm/work_7_20220317/'
adata = dyn.read_h5ad(data_path+"zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")

umap1 = adata.obs['umap_1']
umap2 = adata.obs['umap_2']


fig_path = 'E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/figs_phdThesis_1/vf_construct_1/'

save_fig = True

if save_fig:
    save_show_or_return = 'save'
else:
    save_show_or_return = 'return'


basis = 'umap'

save_kwargs = {"path": fig_path, "prefix": 'cell_wise_vel_ratio_'+basis+'_1', "dpi": None, "ext": 'pdf',
            "transparent": True, "close": True, "verbose": True} 

basis = "umap"
V = adata.obsm[f"velocity_{basis}"]
speed = np.linalg.norm(V, axis=1)

ratio = 0.1 # 只取每个cell_type中，比例为ratio个细胞画cell wise velocity
idx_all = []

for ct, ix in adata.obs.groupby("Cell_type").indices.items():
    ix = np.array(list(ix))
    k_per_type = round(len(ix)*ratio)
    top = ix[np.argsort(speed[ix])[-min(k_per_type, len(ix)) :]]
    idx_all.append(top)

idx = np.unique(np.concatenate(idx_all))

dyn.pl.cell_wise_vectors(
    adata,
    color=['Cell_type'],
    basis=basis,
    cell_inds=idx.tolist(),
    quiver_length=15,
    quiver_size=10,
    pointsize=0.1,
    show_legend='on data',
    save_show_or_return=save_show_or_return,
    save_kwargs=save_kwargs
)


basis = 'pca'
save_kwargs = {"path": fig_path, "prefix": 'cell_wise_vel_ratio_'+basis+'_1', "dpi": None, "ext": 'pdf',
            "transparent": True, "close": True, "verbose": True} 

V = adata.obsm[f"velocity_{basis}"]
speed = np.linalg.norm(V, axis=1)

idx_all = []

for ct, ix in adata.obs.groupby("Cell_type").indices.items():
    ix = np.array(list(ix))
    k_per_type = round(len(ix)*ratio)
    top = ix[np.argsort(speed[ix])[-min(k_per_type, len(ix)) :]]
    idx_all.append(top)

idx = np.unique(np.concatenate(idx_all))

dyn.pl.cell_wise_vectors(
    adata,
    color=['Cell_type'],
    basis=basis,
    cell_inds=idx.tolist(),
    quiver_length=15,
    quiver_size=10,
    pointsize=0.1,
    show_legend='on data',
    save_show_or_return=save_show_or_return,
    save_kwargs=save_kwargs
)


print()