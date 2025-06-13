# 3D Hazard Space
# (Fig.2 C-E)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LogNorm

cmap_name = 'turbo' #'viridis'

Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
# fname = f"{Program_Dir_Path}/LAB_3Dspace.txt"
# outname = f"{Program_Dir_Path}/LAB_3Dspace.png"
# shape_dict = {'UE': 'o', 'WE': 's',
#               'UAF': '*', 'WAF': 'd',
#               'FL': '^'}

fname = f"{Program_Dir_Path}/SRA_3Dspace.txt"
outname = f"{Program_Dir_Path}/SRA_3Dspace.png"
shape_dict = {'HO': 'P', 'HF': '*',
              'MS': 'o', 'TS': 's',
              'AS': 'd', 'CSW': '^',
              'LFW': 'X'}


df = pd.read_csv(fname, sep='\t', index_col=0)
raw_RI = df["RiskIndex"]
min_RI, max_RI = raw_RI.min(), raw_RI.max()
norm = Normalize(vmin=min_RI, vmax=max_RI)
# norm = LogNorm(vmin=min_RI, vmax=max_RI)


fig = plt.figure(figsize=(10, 8), dpi=300)
ax = fig.add_subplot(111, projection='3d')
handle_ls = []

for group_name in shape_dict:
    group_df = df[df["Group"] == group_name]
    x = group_df["ARGs"]
    y = group_df["MGEs"]
    z = group_df["VFs"]
    sizes = np.log10(group_df["CoocurScore"]) * 500
    colors = group_df["RiskIndex"]
    marker = shape_dict[group_name]
    
    scatter = ax.scatter(
        x, y, z,
        c=colors,
        s=sizes,
        alpha=0.7,
        cmap=cmap_name,
        norm=norm,
        marker=marker,
        label=group_name,
        zorder=10,
    )
    
    sub_marker_size = 4
    sub_marker_alpha = 0.7
    sub_marker_color = 'grey'
    scatter_xy = ax.scatter(x, y, 0, s=sub_marker_size, c=sub_marker_color, marker=marker, label=group_name, alpha=sub_marker_alpha, zorder=9)
    scatter_yz = ax.scatter(0, y, z, s=sub_marker_size, c=sub_marker_color, marker=marker, label=group_name, alpha=sub_marker_alpha, zorder=9)
    scatter_xz = ax.scatter(x, 0, z, s=sub_marker_size, c=sub_marker_color, marker=marker, label=group_name, alpha=sub_marker_alpha, zorder=9)
    
    mksize = 14 if marker == '*' else 8
    proxy = Line2D([0], [0], 
           marker=marker, 
           color='w',
           markerfacecolor='k',
           markersize=mksize, 
           label=group_name),
    
    handle_ls.append(proxy[0])


ax.legend(handles=handle_ls, loc='upper right', fontsize=10)

sm = ScalarMappable(norm=norm, cmap=cmap_name)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.05, shrink=0.75)
cbar.set_label('Risk Index', fontsize=14)

ax.set_xlabel('ARG', fontsize=14)
ax.set_ylabel('MGE', fontsize=14)
ax.set_zlabel('VF', fontsize=14)
ax.set_xlim(left=.0)
ax.set_ylim(bottom=.0)
ax.set_zlim(bottom=.0)

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_zscale('log')

# ax.set_title('3D Hazard Space')
ax.set_box_aspect((10, 10, 10))
ax.zaxis.labelpad = 2
ax.view_init(elev=15, azim=50)

plt.tight_layout()
fig.subplots_adjust(left=0.05, right=1.0, bottom=0.00, top=1.0)
plt.savefig(outname, dpi=300)
plt.show()
