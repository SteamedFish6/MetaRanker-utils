# Barplot: High Risk Genes, BPM & Co-ocurrence Frequency (ARGs & MGEs)
# (Fig.5 D)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/HRG_ARG-MGE_category.txt"
outname = f"{Program_Dir_Path}/HRG_ARG-MGE_category.png"

df = pd.read_csv(fname, sep='\t')
data1 = df["BPM"].values
data2 = df["coocur"].values
labels = df['category'].to_list()

fig = plt.figure(figsize=(9, 6.5), dpi=300)
ax = plt.gca()

bar_width = 0.37
length = len(labels)
b1 = ax.bar(x=np.arange(length), height=data1, label='Abundance', color='#17BECF', alpha=0.9, width=bar_width, edgecolor='k', zorder=5)

ax1 = ax.twinx()
b2 = ax1.bar(x=np.arange(length)+bar_width*1, height=data2, label='Co-ocurrence', color='#E377C2', alpha=0.9, width=bar_width, edgecolor='k', zorder=5)

ax.set_ylabel('Abundance (BPM)', fontsize=18)
ax.tick_params(axis='y', labelsize=14)
ax1.set_ylabel('Co-ocurrence Frequency', fontsize=18)
ax1.tick_params(axis='y', labelsize=14)

plt.xlabel('')
ax.set_xticks(ticks=np.arange(length)+bar_width*0.5, labels=labels, fontsize=12, rotation=45)
ax.legend(handles=[b1, b2], fontsize=12, loc='upper right')

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
