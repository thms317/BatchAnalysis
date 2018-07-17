import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters\\clean\\"
save_folder = "C:\\Users\\brouw\\Desktop\\"

df_k = pd.read_csv(folder+"k_pars_all.dat", sep='\t')
# df_G1 = pd.read_csv(folder+"G1_pars_all.dat", sep='\t')
# df_G2 = pd.read_csv(folder+"G2_pars_all.dat", sep='\t')

keys = list(df_k.keys())

boxplot = []
for n,key in enumerate(keys):
    data = df_k[key]
    data = np.array(data[~np.isnan(data)])
    boxplot.append(data)

linewidth = 2
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
bp = ax.boxplot(boxplot)
ax.set_ylim(0.1,3)
ax.semilogy()

for box in bp['boxes']:
    box.set(linewidth=1.5)
## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(linewidth=linewidth)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(linewidth=linewidth)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#b2df8a', linewidth=linewidth)

# ## change the style of fliers and their fill
# for flier in bp['fliers']:
#     flier.set(marker='o', color='#e7298a', alpha=0.5)

ax.set_xticklabels(keys, rotation=45)
ax.set_title("Stiffness (pN/nm)")


fig.savefig(save_folder+"boxplot.png")

plt.close()