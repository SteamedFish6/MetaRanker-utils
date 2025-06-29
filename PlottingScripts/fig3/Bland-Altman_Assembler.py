# Bland-Altman Figure: RI, Megahit vs MetaSpades
# (Fig.3 E)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import pingouin as pg
import matplotlib.pyplot as plt


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
fname = f"{Program_Dir_Path}/Assembler.txt"
outname = f"{Program_Dir_Path}/Assembler_BAplot.png"

g1_name = 'Megahit'
g2_name = 'MetaSpades'

df = pd.read_csv(fname, sep='\t')
group1 = df[g1_name].values
group2 = df[g2_name].values
labels1 = df['SampleName'].to_list()
labels2 = df['Environment'].to_list()
# labels = ["{}, {}".format(lb1, lb2) for (lb1, lb2) in zip(labels1, labels2)]
labels = labels1

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
plt.axhline(diff_mean, color='#D62728', linestyle='--', zorder=5, lw=1.2)
plt.axhline(diff_sd1, color='grey', zorder=5, lw=1)
plt.axhline(diff_sd2, color='grey', zorder=5, lw=1)

xmin, xmax, ymin, ymax = plt.axis()
plt.text(xmin+(xmax-xmin)*0.835, diff_mean+(ymax-ymin)*0.01, 'Mean ({:.2f})'.format(diff_mean), fontsize=14, zorder=8, color='#D62728')
plt.text(xmin+(xmax-xmin)*0.78, diff_sd1+(ymax-ymin)*0.01, '+1.96SD ({:.2f})'.format(diff_sd1), fontsize=14, zorder=8, color='grey')
plt.text(xmin+(xmax-xmin)*0.781, diff_sd2+(ymax-ymin)*0.01, '-1.96SD ({:.2f})'.format(diff_sd2), fontsize=14, zorder=8, color='grey')

plt.text(xmin+(xmax-xmin)*0.04, ymin+(ymax-ymin)*0.8, 'ICC(3, 1) = {:.4f}\np = {:.2e}'.format(icc_value, pvalue), fontsize=16, zorder=9, color='k')

# plt.xlabel(u'$\\mathregular{(RI_{Megahit} + RI_{MetaSpades})/2}$')
# plt.ylabel(u'$\\mathregular{RI_{Megahit} - RI_{MetaSpades}}$')
plt.xlabel('(RI$_{Megahit}$ + RI$_{MetaSpades}$)/2', fontsize=18)
plt.ylabel('RI$_{Megahit}$ - RI$_{MetaSpades}$', fontsize=18)
plt.tick_params(labelsize=14)

plt.tight_layout()
plt.savefig(outname, dpi=300)
plt.show()
