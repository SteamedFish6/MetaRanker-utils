# Regression Analysis: Risk Vector vs RPM, RPKM, BPM (ARGs, MGEs, VFs)
# (Fig.3 A-C)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pingouin as pg


RE_type = "ARGs" # #"MGEs" #"VFs"
Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Abund_validation.txt"
outname = f"{Program_Dir_Path}/Abund_validation_{RE_type}.png"

df = pd.read_csv(fname, sep='\t')
xname_list = ["RPM", "RPKM", "BPM"]
yname_list = [f"{xname} of {RE_type}" for xname in xname_list]
yname_list2 = [f"Adj. {xname} of {RE_type}" for xname in xname_list]

mark_list = ['o', '^', 's']
color_list = ['#FF7F0E', '#2CA02C', '#D62728']
plt.figure(figsize=(9, 6), dpi=300)

text_ls = []
text_ls2 = []
for i in range(len(xname_list)):
    n_samples = df.shape[0]
    x = df[RE_type].values
    y = df[yname_list[i]].values
    y_adj = df[yname_list2[i]].values
    
    coefficients = np.polyfit(x, y, 1)
    slope, intercept = coefficients
    sign = '+' if intercept >= 0 else '-'
    regression_line = np.poly1d(coefficients)
    y_pred = regression_line(x)
    r_squared = r2_score(y, y_pred)
    
    plt.scatter(x, y, color=color_list[i], alpha=0.9, s=36, zorder=10, marker=mark_list[i])
    plt.plot(x, y_pred, color=color_list[i], linewidth=1.8, alpha=0.9, zorder=5, label=yname_list[i])
    
    text1 = f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\nR² = {r_squared:.4f}'
    
    ICC_data = pd.DataFrame({
        'Measurements': np.concatenate([x, y_adj]),
        'Rater': [RE_type for j in range(n_samples)] + [yname_list[i] for j in range(n_samples)],
        'Target': np.arange(n_samples).tolist() * 2
    })
    icc = pg.intraclass_corr(
        data=ICC_data,
        targets='Target',
        raters='Rater',
        ratings='Measurements'
    )
    icc_value = icc[icc['Type'] == 'ICC3']['ICC'].values[0]
    pvalue = icc[icc['Type'] == 'ICC3']['pval'].values[0]
    print(f"ICC(3, 1): {icc_value}, p = {pvalue}")
    text2 = 'ICC(3, 1) = {:.4f}\np = {:.2e}'.format(icc_value, pvalue)
    
    text_ls.append(text1)
    text_ls2.append(text2)

y_sep = 0.1
xmin, xmax, ymin, ymax = plt.axis()
for i in range(len(xname_list)):
    plt.text(xmin+(xmax-xmin)*0.65, ymin+(ymax-ymin)*(0.24-i*y_sep), text_ls[i], fontsize=16, color=color_list[i])
    plt.text(xmin+(xmax-xmin)*0.02, ymin+(ymax-ymin)*(0.68-i*y_sep), text_ls2[i], fontsize=16, color=color_list[i])

plt.legend(fontsize=16, loc="upper left", )
plt.xlabel(f'{RE_type} Value in Risk Vector', fontsize=18)
plt.ylabel(f'{RE_type} Abundance', fontsize=18)

plt.tick_params(labelsize=14)
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
