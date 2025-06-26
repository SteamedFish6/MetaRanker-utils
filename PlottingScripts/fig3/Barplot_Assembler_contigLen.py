# Barplot: contigs length, Megahit vs MetaSpades
# (Fig.3 F)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Assembler_contigLen.txt"
outname = f"{Program_Dir_Path}/Assembler_contigLen.png"

df = pd.read_csv(fname, sep='\t')
data1 = df["Megahit AvgContigLength"].values
data2 = df["MetaSpades AvgContigLength"].values
labels1 = df['SampleName'].to_list()
labels2 = df['Environment'].to_list()
# labels = ["{}\n{}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]
labels = labels1

plt.figure(figsize=(9, 6), dpi=300)
ax = plt.gca()

bar_width = 0.4
length = len(labels)
ax.bar(x=np.arange(length), height=data1, label='Megahit', color='#1F77B4', alpha=0.8, width=bar_width, edgecolor='k', zorder=5)
ax.bar(x=np.arange(length)+bar_width, height=data2, label='MetaSpades', color='#FF7F0E', alpha=0.8, width=bar_width, edgecolor='k', zorder=5)

ax.set_ylabel('Average Length of Contigs (bp)', fontsize=18)
ax.tick_params(axis='y', labelsize=14)

ax.set_xticks(ticks=np.arange(length)+bar_width*0.5, labels=labels, fontsize=11, rotation=90)
ax.legend(fontsize=14, loc="upper right")

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
