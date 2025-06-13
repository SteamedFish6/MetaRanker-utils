# Box-Violin Plot, RI & Adj. Ecological & Adj. Human Health
# (Fig.4 A)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# from scipy import stats


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/MetaComp.txt"
outname = f"{Program_Dir_Path}/MetaComp_box_violin.png"

x_name = "Group"
y_name = "Risk Index"

raw_df = pd.read_csv(fname, sep='\t')
categories = ["MetaRanker", "Adj. Ecological", "Adj. Human Health"]
df = pd.DataFrame()
stdev_dict = {}
for group in categories:
    group_values = raw_df[group].dropna(inplace=False)
    stdev_dict[group] = np.std(group_values)
    temp_df = pd.DataFrame({y_name: group_values, x_name: group})
    df = pd.concat([df, temp_df])

palette = sns.color_palette("pastel", n_colors=len(categories))

fig, ax = plt.subplots(figsize=(9, 6), dpi=300)

violin = sns.violinplot(x=x_name, 
                        y=y_name, 
                        hue=x_name,
                        data=df,
                        inner=None,
                        density_norm='width',
                        width=0.4,
                        linewidth=1,
                        saturation=0.6,
                        palette=palette,
                        dodge=False,
                        zorder=5,
                        ax=ax)

box = sns.boxplot(x=x_name, 
                  y=y_name, 
                  hue=x_name,
                  data=df,
                  width=0.15,
                  linewidth=2,
                  showfliers=False,
                  flierprops={"marker": "x",
                             "markersize": 6,
                             "markerfacecolor": "black",
                             "zorder": 16},
                  boxprops={"facecolor": "grey", "edgecolor": "black", "linewidth": 0.9, "zorder": 12},
                  whiskerprops={"color": "black", "linewidth": 0.9, "zorder": 11},
                  medianprops={"color": "black", "linewidth": 0.9, "zorder": 13},
                  capprops={"color": "black", "linewidth": 0.9, "zorder": 15},
                  showmeans=True,
                  meanprops={
                            #  "marker": 'o',
                             "markersize": 7,
                             "markerfacecolor": "white",
                             "markeredgecolor": "black",
                             "zorder": 14},
                  palette=palette,
                  dodge=False,
                  zorder=10,
                  ax=ax
                  )


for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)

plt.grid(False)
legend = ax.legend(handles=[], frameon=False)

for group in categories:
    x = df[x_name].unique().tolist().index(group)
    y = 500
    plt.text(x, y, "std = {:.2f}".format(stdev_dict[group]), ha='center', va='bottom', color="black", fontsize=14)

plt.ylim(top=550)

plt.xlabel(' ')
plt.ylabel(y_name, fontsize=18)
plt.xticks(ticks=range(len(categories)), labels=categories, fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(outname, dpi=300)
plt.show()

