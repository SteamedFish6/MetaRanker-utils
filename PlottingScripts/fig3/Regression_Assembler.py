# Regression Analysis: RI, Megahit vs MetaSpades
# (Fig.3 D)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Assembler.txt"
outname = f"{Program_Dir_Path}/Assembler_regression.png"

df = pd.read_csv(fname, sep='\t')
x = df['Megahit'].values
y = df['MetaSpades'].values
labels2 = df['SampleName'].to_list()
labels1 = df['Environment'].to_list()
labels = ["{}, {}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]

coefficients = np.polyfit(x, y, 1)
slope, intercept = coefficients
sign = '+' if intercept >= 0 else '-'
regression_line = np.poly1d(coefficients)
y_pred = regression_line(x)
r_squared = r2_score(y, y_pred)

plt.figure(figsize=(9, 8), dpi=300)
plt.scatter(x, y, color='#1F77B4', alpha=0.9, s=36, zorder=10)
plt.plot(x, y_pred, color='#D62728', linewidth=2, alpha=0.9, zorder=5, 
         label=f'Fit: y = {slope:.2f}x + {intercept:.2f}\nR² = {r_squared:.4f}')
plt.xlabel('RI$_{Megahit}$', fontsize=18)
plt.ylabel('RI$_{MetaSpades}$', fontsize=18)
plt.text(8, 280, f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\nR² = {r_squared:.4f}', fontsize=16)

x_bias = [0 for i in range(df.shape[0])]
y_bias = [0 for i in range(df.shape[0])]
# x_bias = [-9, 3.5, -108, 4, 4, 4, 4]
# y_bias = [-9, -3, -3, -3, -4, -5, -7]
for i in range(len(labels)):
    plt.text(x[i]+x_bias[i], y[i]+y_bias[i], labels[i], fontsize=10, zorder=8, color='grey')

plt.tick_params(labelsize=14)
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
