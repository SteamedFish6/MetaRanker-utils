import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

fname = "fig2D-F.txt"
group = 3 #1 #2 #3 #"human feces" #"common soil and water" #"municipal sewage" 
outname = "DataSizeTest_{}_3.png".format(group)
anchor = (0, 1) #(0.2, 0.6)  #(0.98, 0.78)
l_loc = "upper left"
y_lim = 7.75 #140 #500 #140 #7.5

df = pd.read_csv(fname, sep='\t')
df_group = df[df["Group"] == group]

xname_list = df_group["Sample"].unique().tolist()
mark_list = ['s', 'o', '^', 'v']
color_list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728']
plt.figure(figsize=(9, 6), dpi=300)

for i in range(len(xname_list)):
    xname = xname_list[i]
    df_sample = df_group[df_group["Sample"] == xname]
    group_name = df_sample["GroupName"].tolist()[0]
    
    df_sample = df_sample.sort_values(by="DataSize (G)")
    x = df_sample["DataSize (G)"]
    y_all = df_sample[[f"RI1", f"RI2", f"RI3"]]
    y = y_all.sum(axis=1) / 3
    y_error = y_all.std(axis=1, ddof=0)

    x_smooth = np.linspace(x.min(), x.max(), 1000)
    spline = make_interp_spline(x, y, k=3)
    y_smooth = spline(x_smooth)
    
    plt.scatter(x, y, color=color_list[i], alpha=0.8, s=40, zorder=10, marker=mark_list[i], label=f"{xname} ({group_name})")
    plt.errorbar(x, y, yerr=y_error, linewidth=0, alpha=0.8, color=color_list[i], ecolor=color_list[i], capsize=8, capthick=1.5, elinewidth=1.5, zorder=5)
    plt.plot(x_smooth, y_smooth, color=color_list[i], linewidth=1.75, alpha=0.8, zorder=8)

plt.legend(fontsize=14, loc=l_loc, bbox_to_anchor=anchor)
plt.xlabel("DataSize (G)", fontsize=18)
plt.ylabel("RI", fontsize=18)
plt.ylim(top=y_lim)

plt.tick_params(labelsize=14)
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(outname,dpi=300)
plt.show()
