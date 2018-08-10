import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

datafile = "assembled_pars.dat"
folder = "C:\\Users\\tbrouwer\\Desktop\\Chromatin Force Spectroscopy\\"
df = pd.read_csv(folder + datafile, sep='\t')

boxplot, int_keys, sub_keys = [], [], []

NRLs = df.NRL.unique()

for NRL in NRLs:
    stiffness = np.array(df.loc[df["NRL"] == NRL, 'k_pN_nm'])

    boxplot.append(stiffness)
    int_keys.append((int(NRL)))
    sub_keys.append(str(int(NRL))+"x16")

    x = np.random.normal(NRL, 0.05, size=len(stiffness))
    plt.plot(x, stiffness, '.', color="black", zorder=10, alpha=0.2)

plt.boxplot(boxplot, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=45)
plt.title("Stiffness (pN/nm)", fontsize=20)
plt.savefig(folder+"stiffness",dpi=600)
plt.show()

