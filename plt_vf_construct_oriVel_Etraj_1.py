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

adata: Anndata = dyn.read("zebrafish_H44atoh7_1/data/H44DL1208_pres_dynamo_oriVer_processed.h5ad")
# adata: Anndata = dyn.read("zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")



fig_path = 'E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/figs_phdThesis_1/vf_construct_1/'

vfu = adata.uns['VecFld_umap']
Etraj_vfu = vfu['E_traj']

vfpca = adata.uns['VecFld_pca']
Etraj_pca = vfpca['E_traj']


### Etraj_vfu umap plot
# -------------------------
# plt.figure()
# plt.plot(np.arange(1,Etraj_vfu.shape[0]+1),Etraj_vfu,'-',linewidth=2,alpha=0.8)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.xlabel(r'Iterations',fontsize=16)
# plt.ylabel(r'Loss Funtion',fontsize=16)
# # plt.legend(loc='upper right',fontsize=16)
# # plt.savefig(fig_path+'lossFun_umap_oriVer_1.pdf')
# plt.show()

# -------------------------
### Etraj_vfu umap plot, 纵坐标指数坐标
# -------------------------
import matplotlib.ticker as mticker

fig, ax = plt.subplots(figsize=(6.4, 4.8))
ax.plot(np.arange(1, Etraj_vfu.shape[0] + 1), Etraj_vfu-Etraj_vfu.min()+1, '-', linewidth=2, alpha=0.8)

# 关键：对数纵轴
ax.set_yscale('log')

# 让主刻度显示为 10^n（截图风格）
ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10))
ax.minorticks_off()
ax.set_xlabel('Iterations', fontsize=20)
ax.set_ylabel('Loss Function', fontsize=20)
ax.tick_params(axis='both', labelsize=20)
plt.savefig(fig_path+'lossFun_umap_oriVer_2.pdf')
plt.show()


print()