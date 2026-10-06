import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dynamo as dyn
import seaborn as sns
import anndata as Anndata
# type set
from typing import Dict, List, Any, Tuple
from pandas import DataFrame

dyn.dynamo_logger.main_silence()

# filter warnings for cleaner tutorials
import warnings

warnings.filterwarnings('ignore')

adata: Anndata = dyn.read("zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")

# dyn.vf.VectorField(adata, basis='umap', n=100, M=1000,pot_curl_div=True)
dyn.pl.topography(adata, basis='umap', background='white', color=['Cell_type'],streamline_color='black',show_legend='on data', frontier=True,save_show_or_return='return')
# black: absorbing fixed points;
# red: emitting fixed points;
# blue: unstable fixed points.


fig_path = 'E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/figs_phdThesis_1/vf_construct_1/'

### loss function vs iterations
vfu = adata.uns['VecFld_umap']
Etraj_vfu = vfu['E_traj']

vfpca = adata.uns['VecFld_pca']
Etraj_pca = vfpca['E_traj']

### select good fixed points
vfu = adata.uns["VecFld_umap"]
Xss, ftype = adata.uns["VecFld_umap"]["Xss"], adata.uns["VecFld_umap"]["ftype"]
good_fixed_points = [60, 27, 14, 24, 22, 25, 36, 12, 48, 33, 37, 40,3] 
# good_fixed_points = np.where(ftype==0)[0]
# good_fixed_points = np.where(ftype==1)[0]
adata.uns["VecFld_umap"]["Xss"] = Xss[good_fixed_points]
adata.uns["VecFld_umap"]["ftype"] = ftype[good_fixed_points]
# s_id = 1
# save_tmp= {"path": None, "prefix": 'zebrafish_H44atoh7_1/Figures_H44_DL_1208_pres/umap_adata/vectorfield_plot_tol_redu_stable_'+str(s_id), "dpi": None, "ext": 'pdf',
#             "transparent": True, "close": True, "verbose": True}
# dyn.pl.topography(adata, basis='umap', background='white', color=['Cell_type'], streamline_color='black',show_legend='on data', frontier=True, save_show_or_return='save',save_kwargs=save_tmp)

dyn.pl.topography(adata, basis='umap', background='white', color=['Cell_type'], streamline_color='black',show_legend='on data', frontier=True)


### Etraj_vfu umap plot
# -------------------------

# plt.figure()
# plt.plot(np.arange(1,Etraj_vfu.shape[0]+1),Etraj_vfu,'-',linewidth=2,alpha=0.8)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.xlabel(r'Iterations',fontsize=20)
# plt.ylabel(r'Loss Funtion',fontsize=20)
# # plt.legend(loc='upper right',fontsize=20)
# plt.savefig(fig_path+'lossFun_umap_newVer_2.pdf')
# plt.show()

# -------------------------
### Etraj_vfu umap plot, 纵坐标指数坐标
# -------------------------
import matplotlib.ticker as mticker

fig, ax = plt.subplots(figsize=(5.8, 4.8))
ax.plot(np.arange(1, Etraj_vfu.shape[0] + 1), Etraj_vfu-Etraj_vfu.min()+1, '-', linewidth=2, alpha=0.8)

# 关键：对数纵轴
ax.set_yscale('log')

# 让主刻度显示为 10^n（截图风格）
ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10))

# # （可选）加一些次刻度：2,3,...,9 × 10^n
# ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
# ax.yaxis.set_minor_formatter(mticker.NullFormatter())
ax.minorticks_off()
ax.set_xlabel('Iterations', fontsize=20)
ax.set_ylabel('Loss Function', fontsize=20)
ax.tick_params(axis='both', labelsize=20)
plt.savefig(fig_path+'lossFun_umap_newVer_2.pdf')
plt.show()


print()