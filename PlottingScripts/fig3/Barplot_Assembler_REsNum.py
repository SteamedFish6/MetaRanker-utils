# Barplot: Detected REs number, Megahit vs MetaSpades
# (Fig.3 H)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Assembler_contigLen.txt"
outname = f"{Program_Dir_Path}/Assembler_contigLen.png"

df = pd.read_csv(fname, sep='\t')
data1 = df["Megahit nREs"].values
data2 = df["MetaSpades nREs"].values
data1_a = df["Megahit nARGs"].values
data2_a = df["MetaSpades nARGs"].values
data1_m = df["Megahit nMGEs"].values
data2_m = df["MetaSpades nMGEs"].values
data1_v = df["Megahit nVFs"].values
data2_v = df["MetaSpades nVFs"].values
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
ax.bar(x=np.arange(length), height=data1_a, label='ARGs', hatch='//', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length)+bar_width, height=data2_a, hatch='//', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length), height=data1_m, bottom=data1_a, label='MGEs', hatch='\\\\', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length)+bar_width, height=data2_m, bottom=data2_a, hatch='\\\\', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length), height=data1_v, bottom=data1_m+data1_a, label='VFs', hatch='xx', fill=False, width=bar_width, edgecolor='k', zorder=8)
ax.bar(x=np.arange(length)+bar_width, height=data2_v, bottom=data2_m+data2_a, hatch='xx', fill=False, width=bar_width, edgecolor='k', zorder=8)

ax.set_ylabel('Number of REs', fontsize=18)
ax.tick_params(axis='y', labelsize=14)

ax.set_xticks(ticks=np.arange(length)+bar_width*0.5, labels=labels, fontsize=11, rotation=90)
ax.legend(fontsize=14, loc="upper right")

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
