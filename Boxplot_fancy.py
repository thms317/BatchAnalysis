import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters\\clean\\"
save_folder = "C:\\Users\\tbrouwer\\Desktop\\Boxplots\\"

df_k = pd.read_csv(folder+"k_pars_all.dat", sep='\t')
df_G1 = pd.read_csv(folder+"G1_pars_all.dat", sep='\t')
df_G2 = pd.read_csv(folder+"G2_pars_all.dat", sep='\t')


NRLs = [3, 4, 6, 8, 10, 11, 12]
NRLs = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

# k

plt.figure(figsize=(12,6))
data = []
keys = list(df_k.keys())
for n, key in enumerate(keys):
    sub_data = df_k[key]
    sub_data = np.array(sub_data[~np.isnan(sub_data)])
    data.append(sub_data)

boxplot, sub_keys, int_keys, length = [], [], [], []
k_tot = np.array([])

for n, i in enumerate(NRLs):
    boxplot.append(data[i])
    sub_keys.append(keys[i])
    NRL = int(keys[i][:3])
    int_keys.append(NRL)
    k_tot = np.concatenate((k_tot, data[i]))
    length.append(len(data[i]))

    # statistic, norm_p_val = stats.normaltest(data[i])
    # print(str(keys[i])+" - normaltest value  =", norm_p_val)

    x = np.random.normal(NRL, 0.05, size=len(data[i]))
    plt.plot(x, data[i], '.', color="black", zorder=10, alpha=0.2)

plt.boxplot(boxplot, showfliers=False, positions=int_keys)
plt.xticks(int_keys,sub_keys, rotation=45)
plt.title("Stiffness (pN/nm)", fontsize=20)
plt.savefig(save_folder+"k_boxplot_fancy",dpi=600)
# plt.show()
plt.close()

# G1

plt.figure(figsize=(12,6))
data = []
keys = list(df_G1.keys())
for n, key in enumerate(keys):
    sub_data = df_G1[key]
    sub_data = np.array(sub_data[~np.isnan(sub_data)])
    data.append(sub_data)

boxplot, sub_keys, int_keys, length = [], [], [], []
G1_tot = np.array([])

for n, i in enumerate(NRLs):
    boxplot.append(data[i])
    sub_keys.append(keys[i])
    NRL = int(keys[i][:3])
    int_keys.append(NRL)
    G1_tot = np.concatenate((G1_tot, data[i]))
    length.append(len(data[i]))

    # statistic, norm_p_val = stats.normaltest(data[i])
    # print(str(keys[i])+" - normaltest value  =", norm_p_val)

    x = np.random.normal(NRL, 0.05, size=len(data[i]))
    plt.plot(x, data[i], '.', color="black", zorder=10, alpha=0.2)

plt.boxplot(boxplot, showfliers=False, positions=int_keys)
plt.xticks(int_keys,sub_keys, rotation=45)
plt.title("G1 (kT)", fontsize=20)
plt.savefig(save_folder+"G1_boxplot_fancy",dpi=600)
# plt.show()
plt.close()

# G2

plt.figure(figsize=(12,6))
data = []
keys = list(df_G2.keys())
for n, key in enumerate(keys):
    sub_data = df_G2[key]
    sub_data = np.array(sub_data[~np.isnan(sub_data)])
    data.append(sub_data)

boxplot, sub_keys, int_keys, length = [], [], [], []
G2_tot = np.array([])

for n, i in enumerate(NRLs):
    boxplot.append(data[i])
    sub_keys.append(keys[i])
    NRL = int(keys[i][:3])
    int_keys.append(NRL)
    G2_tot = np.concatenate((G2_tot, data[i]))
    length.append(len(data[i]))

    # statistic, norm_p_val = stats.normaltest(data[i])
    # print(str(keys[i])+" - normaltest value  =", norm_p_val)

    x = np.random.normal(NRL, 0.05, size=len(data[i]))
    plt.plot(x, data[i], '.', color="black", zorder=10, alpha=0.2)

plt.boxplot(boxplot, showfliers=False, positions=int_keys)
plt.xticks(int_keys,sub_keys, rotation=45)
plt.title("G2 (kT)", fontsize=20)
plt.savefig(save_folder+"G2_boxplot_fancy",dpi=600)
# plt.show()
plt.close()