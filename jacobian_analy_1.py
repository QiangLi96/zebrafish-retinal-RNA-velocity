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
adata: Anndata = dyn.read("E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/data/H44DL1208_pres_processed_data.h5ad")


fig_path = 'E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/phdThesis_1/jacibian_1/'

selected_genes = ['ptf1a','tfap2a']
# selected_genes = ['atoh7','ptf1a']

# 裁剪 colormap（推荐）
import matplotlib.colors as colors
base = plt.get_cmap("YlGn")
cmap = colors.LinearSegmentedColormap.from_list(
    "YlGn_cut", base(np.linspace(0, 0.8, 256))   # 可改小/改大
)
# 用你的绘图对象里把 cmap="YlGn" 换成 cmap=cmap
# 例如：
# surf = ax.plot_surface(..., cmap=cmap, ...)
# 或 imshow/pcolormesh 也一样

fig_path = 'E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/figs_phdThesis_1/jacobian_1/'
fig_name = 'genes_expre_ptf1a_tfap2a_1'
# fig_name = 'genes_expre_atoh7_ptf1a_1'

save_kwargs= {"path": fig_path, "prefix": fig_name, "dpi": None, "ext": 'pdf', "transparent":
            True, "close": True, "verbose": True}

dyn.pl.umap(adata,cmap=cmap, color=selected_genes, layer='M_s', frontier=True,save_show_or_return='save',save_kwargs=save_kwargs)


dyn.vf.jacobian(adata, regulators=selected_genes, effectors=selected_genes)
dyn.pl.jacobian(adata,regulators=selected_genes,effectors=selected_genes, basis="umap")



print()


# plt.figure()
# plt.plot(np.arange(1,Etraj_vfu.shape[0]+1),Etraj_vfu,'-o',linewidth=2,alpha=0.8)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.xlabel(r'Iterations',fontsize=20)
# plt.ylabel(r'Loss Funtion',fontsize=20)
# # plt.legend(loc='upper right',fontsize=20)
# # plt.savefig(fig_path+'lossFun_umap_newVer_1.pdf')
# plt.show()