import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

datafile = "assembled_pars.dat"
folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Analysis\\"
df = pd.read_csv(folder + datafile, sep='\t')

plt.rcParams.update({'font.size': 20})  # legend + title size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

boxplot_k, boxplot_G1, boxplot_G2, int_keys, sub_keys = [], [], [], [], []

NRLs = df.NRL.unique()

for NRL in NRLs:

    stiffness = np.array(df.loc[df["NRL"] == NRL, 'k_pN_nm'])
    G1 = np.array(df.loc[df["NRL"] == NRL, 'G1_kT'])
    G2 = np.array(df.loc[df["NRL"] == NRL, 'G2_kT'])

    boxplot_k.append(stiffness)
    boxplot_G1.append(G1)
    boxplot_G2.append(G2)

    int_keys.append((int(NRL)))
    sub_keys.append(str(int(NRL)))

    x = np.random.normal(NRL, 0.05, size=len(stiffness))
    y = np.random.normal(NRL, 0.05, size=len(G1))
    z = np.random.normal(NRL, 0.05, size=len(G2))

    plt.figure(0, figsize=(12,7))
    plt.plot(x, stiffness, '.', color="black", zorder=10, alpha=0.2)
    plt.figure(1, figsize=(12,7))
    plt.plot(y, G1, '.', color="black", zorder=10, alpha=0.2)
    plt.figure(2, figsize=(12,7))
    plt.plot(z, G2, '.', color="black", zorder=10, alpha=0.2)

plt.figure(0)
plt.boxplot(boxplot_k, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Stiffness (pN/nm)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Stiffness (kN/nm)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
plt.savefig(folder+"stiffness_boxplot",dpi=600)

plt.figure(1)
plt.boxplot(boxplot_G1, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Stacking Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Stacking Energy (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
plt.savefig(folder+"G1_boxplot",dpi=600)

plt.figure(2)
plt.boxplot(boxplot_G2, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Partial Unwrapping Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Partial Unwrapping Energy (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
plt.savefig(folder+"G2_boxplot",dpi=600)

plt.close()

