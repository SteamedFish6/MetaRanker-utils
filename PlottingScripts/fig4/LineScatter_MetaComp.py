# Line-Scatter Plot, RI & Adj. Ecological & Adj. Human Health
# (Fig.4 E)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from scipy import stats


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/MetaComp.txt"
outname = f"{Program_Dir_Path}/MetaComp_line_scatter.png"

x_name = "Group"
y_name = "Risk Index"

df = pd.read_csv(fname, sep='\t')
categories = ["MetaRanker", "Ecological", "Human Health", "Adj. Ecological", "Adj. Human Health"]
color_ls = ['#1F77B4', '#2CA02C', '#D62728', '#2CA02C', '#D62728']
linestyle_ls = ['-', '--', '--', '-', '-']
marker_ls = ['s', '^', '^', 's', 's']
# label_ls = ["Ranker", "Adj. Ecological", "Adj. Human Health", None, None]

x_labels1 = df['Sample'].to_list()
x_labels2 = df['Environment'].to_list()
x_labels = ["{}\n{}".format(lb1, lb2) for (lb1, lb2) in zip(x_labels1, x_labels2)]

x = np.arange(len(x_labels))
plt.figure(figsize=(9, 6), dpi=300)

for i in range(len(categories)):
    plt.scatter(x, df[categories[i]].values, marker=marker_ls[i], color='grey', s=20, zorder=10)
    plt.plot(x, df[categories[i]].values, color=color_ls[i], linestyle=linestyle_ls[i], label=categories[i], linewidth=1.5, zorder=5)

plt.xlabel('')
plt.ylabel(y_name, fontsize=18)
plt.xticks(ticks=range(len(x_labels)), labels=x_labels, rotation=45, fontsize=12)
plt.yticks(fontsize=14)

plt.legend(fontsize=12, loc="upper right")

plt.tight_layout()
plt.savefig(outname, dpi=300)
plt.show()
