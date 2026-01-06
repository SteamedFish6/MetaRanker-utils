import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

Score_Type = "ERS" #"HHRS" "RI"
fname = "fig2.tsv"
outname = "Benchmark" + f"_{Score_Type}.png"

df = pd.read_csv(fname, sep='\t')

xname_list = ["SRR22519794 + QLFW", "HRG + LRG", "SRR22519794 + LRG", "HRG + QLFW"]
mark_list = ['s', 'o', '^', 'v']
color_list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728']
plt.figure(figsize=(9, 6), dpi=300)

text_ls = []
for i in range(len(xname_list)):
    group_name = xname_list[i]
    df_group = df[df["Group"] == group_name]
    x = df_group["PollutionRate"].values
    y_all = df_group[[f"{Score_Type}1", f"{Score_Type}2", f"{Score_Type}3"]]
    y = y_all.sum(axis=1) / 3
    y_error = y_all.std(axis=1, ddof=0)
    y = y.values
    y_error = y_error.values
    
    coefficients = np.polyfit(x, y, 1)
    slope, intercept = coefficients
    sign = '+' if intercept >= 0 else '-'
    regression_line = np.poly1d(coefficients)
    y_pred = regression_line(x)
    
    r_squared = r2_score(y, y_pred)
    
    plt.scatter(x, y, color=color_list[i], alpha=0.8, s=40, zorder=10, marker=mark_list[i])
    plt.errorbar(x, y, yerr=y_error, linewidth=1, linestyle='--', alpha=0.8, color=color_list[i], ecolor=color_list[i], capsize=8, capthick=1.5, elinewidth=1.5, zorder=5)
    plt.plot(x, y_pred, color=color_list[i], linewidth=2, alpha=0.8, zorder=8, label=group_name)
    
    text1 = f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\nR² = {r_squared:.4f}'
    text_ls.append(text1)
    

y_sep = 0.1

plt.ylim(top=155)
xmin, xmax, ymin, ymax = plt.axis()
for i in range(len(xname_list[0:2])):
    plt.text(xmin+(xmax-xmin)*0.4, ymin+(ymax-ymin)*(0.87-i*y_sep), text_ls[i], fontsize=16, color=color_list[i])
xmin, xmax, ymin, ymax = plt.axis()
for i in range(len(xname_list[2:4])):
    plt.text(xmin+(xmax-xmin)*0.7, ymin+(ymax-ymin)*(0.87-i*y_sep), text_ls[i+2], fontsize=16, color=color_list[i+2])

plt.legend(fontsize=14, loc="upper left", )
plt.xlabel('Pollution Level', fontsize=18)
plt.ylabel(Score_Type, fontsize=18)

plt.tick_params(labelsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()

