# Significance Test
# (Fig.2 A, B)
# Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
# fname = f"{Program_Dir_Path}/PlottingScripts/fig2/LAB_data.txt"
# outname = f"{Program_Dir_Path}/PlottingScripts/fig2/LAB_figure_significance.png"

fname = f"{Program_Dir_Path}/PlottingScripts/fig2/SRA_data.txt"
outname = f"{Program_Dir_Path}/PlottingScripts/fig2/SRA_figure_significance.png"


x_name = "Group"
y_name = "Risk Index"

raw_df = pd.read_csv(fname, sep='\t')
categories = list(raw_df.columns)
df = pd.DataFrame()
for group in categories:
    group_values = raw_df[group].dropna(inplace=False)
    temp_df = pd.DataFrame({y_name: group_values, x_name: group})
    df = pd.concat([df, temp_df])

# comparisons = [
#                (categories[0], categories[1]),
#                (categories[2], categories[3]),
#                (categories[3], categories[4]),
#                (categories[1], categories[4]),
#                ]
comparisons = [
               (categories[0], categories[1]),
               (categories[2], categories[3]),
               (categories[0], categories[3]),
               (categories[3], categories[4]),
               (categories[4], categories[5]),
               (categories[5], categories[6]),
               ]

sns.set_theme(style="white", palette="pastel")
palette = sns.color_palette("pastel", n_colors=len(categories))
plt.figure(figsize=(9, 6), dpi=300)

box = sns.boxplot(x=x_name, 
                  y=y_name, 
                  hue=x_name,
                  data=df,
                  width=0.4,
                  linewidth=2,
                  showfliers=False,
                  flierprops={"marker": "x",
                             "markersize": 6,
                             "markerfacecolor": "black",
                             "zorder": 16},
                  boxprops={"edgecolor": "black", "linewidth": 1, "zorder": 12},
                  whiskerprops={"color": "black", "linewidth": 1, "zorder": 11},
                  medianprops={"color": "black", "linewidth": 1, "zorder": 13},
                  capprops={"color": "black", "linewidth": 1, "zorder": 15},
                  showmeans=True,
                  meanprops={
                             "marker": 'x',
                             "markersize": 7,
                             "markerfacecolor": "white",
                             "markeredgecolor": "black",
                             "zorder": 14},
                  palette=palette,
                  dodge=False,
                  zorder=10,
                #   ax=ax,
                  )

# ========== Significance Test ==========
def get_pvalue_stars(p):
    """将p值转换为星号标注"""
    if p < 0.0001: return '****'
    elif p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'

# Shapiro-Wilk
_, p_normal = stats.shapiro(raw_df[categories])

# Levene's Test
groups_tuple = (raw_df[group].dropna(inplace=False) for group in categories)
_, p_var = stats.levene(*groups_tuple)

# Select Method
if p_normal > 0.05 and p_var > 0.05:
    use_method = "parametric"
else:
    use_method = "non-parametric"
print("{}: {}, {}".format(use_method, p_normal, p_var))

def compare_test(group1, group2, method="parametric"):
    if method == "parametric":
        t_stat, p_value = stats.ttest_ind(group1, group2)
    else:
        t_stat, p_value = stats.mannwhitneyu(group1, group2)
    return p_value


y_max = df[y_name].max() * 0.45
default_offset = y_max * 0.05
vertical_offset = default_offset
sub_vertical_offset = y_max * 0.05 * 0.6


for (grp1, grp2) in comparisons:
    
    data1 = df[df[x_name] == grp1][y_name]
    data2 = df[df[x_name] == grp2][y_name]
    
    x1 = df[x_name].unique().tolist().index(grp1)
    x2 = df[x_name].unique().tolist().index(grp2)
    
    p = compare_test(data1, data2, use_method)
    print(p)
    
    current_y = y_max + vertical_offset
    plt.plot([x1, x1, x2, x2], 
             [current_y, current_y+sub_vertical_offset, current_y+sub_vertical_offset, current_y], 
             lw=1.2, color="black")
    
    plt.text((x1+x2)*0.5, current_y+sub_vertical_offset, 
             get_pvalue_stars(p), 
             ha='center', va='bottom', color="black")
    
    vertical_offset += 2*default_offset

plt.ylim(top=y_max + vertical_offset + sub_vertical_offset)
legend = box.legend(handles=[], frameon=False)
plt.xlabel('')
plt.ylabel(y_name, fontsize=18)

plt.xticks(ticks=range(len(categories)), labels=categories, fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(outname, dpi=300)
plt.show()
