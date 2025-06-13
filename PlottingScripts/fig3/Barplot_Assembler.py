# Barplot: RI, Megahit vs MetaSpades
# (Fig.3 F)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Assembler_contigInfo.txt"
outname = f"{Program_Dir_Path}/Assembler_contigNumLength.png"

df = pd.read_csv(fname, sep='\t')
data1 = df["Megahit ContigsNum"].values
data2 = df["MetaSpades ContigsNum"].values
data1_n = df["Megahit N50"].values
data2_n = df["MetaSpades N50"].values
data1_l = df["Megahit AvgContigLength"].values
data2_l = df["MetaSpades AvgContigLength"].values
labels1 = df['SampleName'].to_list()
labels2 = df['Environment'].to_list()
labels = ["{}\n{}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]

plt.figure(figsize=(10, 8), dpi=300)
ax = plt.gca()

bar_width = 0.2
length = len(labels)
ax.bar(x=np.arange(length), height=data1, label='Megahit', color='#1F77B4', alpha=0.9, width=bar_width, edgecolor='k', zorder=5)
ax.bar(x=np.arange(length)+bar_width, height=data2, label='MetaSpades', color='#FF7F0E', alpha=0.9, width=bar_width, edgecolor='k', zorder=5)
ax.bar(x=np.arange(length), height=data1_n, label='N50', hatch='//', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length)+bar_width, height=data2_n, hatch='//', fill=False, width=bar_width, edgecolor='k', zorder=8)

ax1 = ax.twinx()
ax1.bar(x=np.arange(length)+bar_width*2, height=data1_l, label='Megahit', color='#AEC7E8', alpha=1, width=bar_width, edgecolor='k', zorder=5)
ax1.bar(x=np.arange(length)+bar_width*3, height=data2_l, label='MetaSpades', color='#FFBB78', alpha=1, width=bar_width, edgecolor='k', zorder=5)


ax.set_ylabel('Number of Contigs', fontsize=24)
ax.tick_params(axis='y', labelsize=14)
ax.set_ylim([0, 500000])
ax1.set_ylabel('Average Length of Contigs', fontsize=24)
ax1.tick_params(axis='y', labelsize=14)
ax1.set_ylim([500, 3000])

plt.xlabel('')
ax.set_xticks(ticks=np.arange(length)+bar_width*1.5, labels=labels, fontsize=10.5, rotation=30)
# plt.yticks(fontsize=18)
ax.legend(fontsize=16, loc="upper left")
ax1.legend(fontsize=16, loc="upper right")

for spine in ax.spines.values():
    # spine.set_color('black')
    spine.set_linewidth(1.2)
plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
