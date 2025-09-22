# Bland-Altman Figure: RI & Adj. Ecological & Adj. Human Health
# (Fig.4 C, D)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
# from scipy import stats
import pingouin as pg
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/MetaComp.txt"

# outname = f"{Program_Dir_Path}/MetaComp_BAplot_Eco.png"
# g1_name = "MetaRanker"
# g2_name = "Adj. Ecological"
# mid_color = '#2CA02C'

outname = f"{Program_Dir_Path}/MetaComp_BAplot_HH.png"
g1_name = "MetaRanker"
g2_name = "Adj. Human Health"
mid_color = '#D62728'

df = pd.read_csv(fname, sep='\t')
group1 = df[g1_name].values
group2 = df[g2_name].values
labels1 = df['Sample'].to_list()
labels2 = df['Environment'].to_list()
labels = ["{}, {}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]

n_samples = df.shape[0]
data = pd.DataFrame({
    'Measurements': np.concatenate([group1, group2]),
    'Rater': [g1_name for i in range(n_samples)] + [g2_name for i in range(n_samples)],
    'Target': np.arange(n_samples).tolist() * 2
})

icc = pg.intraclass_corr(
    data=data,
    targets='Target',
    raters='Rater',
    ratings='Measurements'
)

icc_value = icc[icc['Type'] == 'ICC3']['ICC'].values[0]
pvalue = icc[icc['Type'] == 'ICC3']['pval'].values[0]
print(f"ICC(3, 1): {icc_value}")

plt.figure(figsize=(9, 6), dpi=300)
mean = (group1 + group2) / 2
diff = group1 - group2
plt.scatter(mean, diff, alpha=0.9, c='#1F77B4', s=36, zorder=10)
diff_mean = np.mean(diff)
diff_sd1 = np.mean(diff) + 1.96 * np.std(diff)
diff_sd2 = np.mean(diff) - 1.96 * np.std(diff)
plt.axhline(diff_mean, color=mid_color, linestyle='--', zorder=5, lw=1)
plt.axhline(diff_sd1, color='grey', zorder=5, lw=1)
plt.axhline(diff_sd2, color='grey', zorder=5, lw=1)

# for i in range(len(diff)):
#     if diff[i] > diff_sd1 or diff[i] < diff_sd2:
#         print(labels[i])
#         # plt.text(mean[i]-61.5, diff[i]-1.85, labels[i], fontsize=8, zorder=7, color='grey')
#         plt.text(mean[i]-67, diff[i]-1.65, labels[i], fontsize=8, zorder=7, color='grey')

xmin, xmax, ymin, ymax = plt.axis()
plt.text(xmin+(xmax-xmin)*0.835, diff_mean+(ymax-ymin)*0.01, 'Mean ({:.2f})'.format(diff_mean), fontsize=14, zorder=8, color=mid_color)
plt.text(xmin+(xmax-xmin)*0.76, diff_sd1+(ymax-ymin)*0.01, '+1.96SD ({:.2f})'.format(diff_sd1), fontsize=14, zorder=8, color='grey') #*0.69, *-0.05
plt.text(xmin+(xmax-xmin)*0.761, diff_sd2+(ymax-ymin)*0.01, '-1.96SD ({:.2f})'.format(diff_sd2), fontsize=14, zorder=8, color='grey') #*0.69

# plt.text(xmin+(xmax-xmin)*0.838, diff_mean+(ymax-ymin)*0.01, 'Mean (0.00)', fontsize=14, zorder=8, color=mid_color)
# plt.text(xmin+(xmax-xmin)*0.76, diff_sd1+(ymax-ymin)*0.01, '+1.96SD ({:.2f})'.format(diff_sd1), fontsize=14, zorder=8, color='grey') #*0.69, *-0.05
# plt.text(xmin+(xmax-xmin)*0.761, diff_sd2+(ymax-ymin)*0.01, '-1.96SD ({:.2f})'.format(diff_sd2), fontsize=14, zorder=8, color='grey') #*0.69

plt.text(xmin+(xmax-xmin)*0.03, ymin+(ymax-ymin)*0.81, 'ICC(3, 1) = {:.4f}\np = {:.2e}'.format(icc_value, pvalue), fontsize=16, zorder=8, color='k')

plt.xlabel(f'(Risk Index + {g2_name})/2', fontsize=18)
plt.ylabel(f'Risk Index - {g2_name}', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(outname, dpi=300)
plt.show()
