import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data_path='E:/python_work_PyCharm/work_7_20220317/zebrafish_H44atoh7_1/figs_phdThesis_1/lap_1/'

# # RGCs parameters
# cell_type = ['RGCs','RGCs']
# start_coord = np.array([[-1.50,2.7],[-1.50,2.7]])
# attract_list=[14,24]
# startCellType=['RGCs','RGCs']

### plot multi LAP in each cell type, to compare different attractors
# # GABAergic ACs parameters
# cell_type = ['GABAergic ACs']
# start_coord = np.array([[-1.35,0.45]])
# startCellType=['GABAergic ACs']
# attract_list=[22]


# # BCs parameters
cell_type = ['BCs','BCs','BCs','BCs']
start_coord = np.array([[2,0.4],[2,0.4],[2,0.4],[2,0.4]])
attract_list=[33,37,40,3]
startCellType=['BCs','BCs','BCs','BCs']

fontsize=14

save_fig = True


for i in range(len(cell_type)):
    kinetics_heatmap = pd.read_csv(data_path+cell_type[i]+'/'+str(attract_list[i])+'_kinetics_heatmap_1.csv',index_col=0)

    
    plt.figure(figsize=(7.0,4.2))
    ax = sns.heatmap(kinetics_heatmap, cmap="viridis",cbar=True)

    # ticks 字体
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    # 纵坐标取消 ticks（含刻度线和刻度标签）
    ax.set_xticks([])
    # colorbar 字体大小
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=fontsize)       # colorbar刻度
    cbar.set_label("normalized expression", fontsize=fontsize)
    plt.tight_layout()
    # # 固定边距（两张图用同一组参数，热图本体大小就一致）
    plt.subplots_adjust(left=0.28, right=0.92, bottom=0.08, top=0.98)  # 关键：固定 left/right
    # plt.subplots_adjust(left=0.08, right=0.92, bottom=0.35, top=0.98)

    if save_fig:
        plt.savefig(data_path+cell_type[i]+'/'+str(attract_list[i])+'_kinetics_heatmap_2.pdf')

plt.show()

print()