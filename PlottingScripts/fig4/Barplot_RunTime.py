# Barplot: Run Time, MeteRanker vs MetaCompare2.0
# (Fig.4 F)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/RunTime.txt"
outname = f"{Program_Dir_Path}/RunTime.png"

df = pd.read_csv(fname, sep='\t')
data1 = df['MetaCompare 2.0'].values
data2 = df['MetaRanker'].values

x_labels1 = df['Sample'].to_list()
x_labels2 = df['Environment'].to_list()
labels = ["{}\n{}".format(lb1, lb2) for (lb1, lb2) in zip(x_labels1, x_labels2)]

plt.figure(figsize=(9, 6.5), dpi=300)
bar_width = 0.4
length = len(labels)
plt.bar(x=np.arange(length), height=data1, label='MetaCompare 2.0', color='steelblue', alpha=0.9, width=bar_width, edgecolor='k')
plt.bar(x=np.arange(length)+bar_width, height=data2, label='MetaRanker', color='indianred', alpha=0.9, width=bar_width, edgecolor='k')

plt.xlabel('')
plt.ylabel('Run Time (s)', fontsize=18)
plt.xticks(ticks=np.arange(length)+bar_width/2, labels=labels, rotation=90, fontsize=11)
plt.yticks(fontsize=14)

plt.legend(fontsize=14, loc="upper left")
plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
