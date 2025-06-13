# Regression Analysis: RI vs Adj. Ecological & Adj. Human Health
# (Fig.4 B)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/MetaComp.txt"
outname = f"{Program_Dir_Path}/MetaComp_regressionplot.png"

df = pd.read_csv(fname, sep='\t')
df = df.sort_values(by='MetaRanker', inplace=True) # sort here
x = df['MetaRanker'].values
y2 = df['Adj. Human Health'].values
y1 = df['Adj. Ecological'].values
labels2 = df['Sample'].to_list()
labels1 = df['Environment'].to_list()
labels = ["{}, {}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]

plt.figure(figsize=(9, 6), dpi=300)
coefficients = np.polyfit(x, y1, 1)
slope, intercept = coefficients
sign = '+' if intercept >= 0 else '-'
regression_line = np.poly1d(coefficients)
y_pred = regression_line(x)
r_squared = r2_score(y1, y_pred)

plt.scatter(x, y1, color='#2CA02C', alpha=0.9, s=30, zorder=10, marker='^')
plt.scatter(x, y2, color='#D62728', alpha=0.9, s=30, zorder=10, marker='o')
plt.plot(x, y_pred, color='#2CA02C', linewidth=1.5, alpha=0.9, zorder=5, label='RI vs Adj. Ecological')
plt.text(240, 40, f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\nR² = {r_squared:.4f}', fontsize=14, color='#2CA02C')

coefficients = np.polyfit(x, y2, 1)
slope, intercept = coefficients
sign = '+' if intercept >= 0 else '-'
regression_line = np.poly1d(coefficients)
y_pred = regression_line(x)
r_squared = r2_score(y2, y_pred)

plt.plot(x, y_pred, color='#D62728', linewidth=1.5, alpha=0.9, zorder=6, label='RI vs Adj. Human Health')
plt.text(240, 10, f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\nR² = {r_squared:.4f}', fontsize=14, color='#D62728')

plt.legend(fontsize=12, loc="upper left")
plt.xlabel('Risk Index', fontsize=18)
plt.ylabel('MetaCompare Adj. Score', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

x_bias = [0 for i in range(df.shape[0])]
y_bias = [0 for i in range(df.shape[0])]
# x_bias = [12, 1, 25, 14,
#           12, 12, -42, -45,
#           -60, -38, -33, -63]
# y_bias = [-7.3, -6.8, -5.8, -1.2,
#           5.5, 1, 10, 4, 
#           4, 5, 4, -9]
for i in range(len(labels)):
    plt.text(x[i]+x_bias[i], y2[i]+y_bias[i], labels[i], fontsize=8, zorder=8, color='grey')

for j in (0, 2, 3, 4, 5):
    plt.plot([x[j], x[j]+x_bias[j]-0.4], [y2[j], y2[j]+y_bias[j]+1.8], lw=0.4, color='k', alpha=1.0)
for j in (6,):
    plt.plot([x[j], x[j]+x_bias[j]+30], [y2[j], y2[j]+y_bias[j]-1.2], lw=0.4, color='k', alpha=1.0)

plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
