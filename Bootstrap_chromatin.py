import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Analysis\\"
datafile = "assembled_pars.dat"

df = pd.read_csv(folder + datafile, sep='\t')

plt.rcParams.update({'font.size': 20})  # legend + title size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

G2_arr, int_keys, sub_keys, length = [], [], [], []

NRLs = df.NRL.unique()

for NRL in NRLs:

    G2 = np.array(df.loc[df["NRL"] == NRL, 'G2_kT'])
    length.append(len(G2))

    G2_arr.append(G2)

    int_keys.append((int(NRL)))
    sub_keys.append(str(int(NRL))+"x16")

# ANOVA
f_val, p_val = stats.f_oneway(G2_arr[0],G2_arr[1],G2_arr[2],G2_arr[3],G2_arr[4],G2_arr[5],G2_arr[6],G2_arr[7],G2_arr[8],G2_arr[9],G2_arr[10])
print("One-way ANOVA P =", p_val)

# remove lowest number of measurements (or not)
# length.remove(min(length))

n_min = min(length)

def bootstrap_resample(X, n):
    X = np.array(X)
    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    X_resample = np.array(X[resample_i])
    return X_resample

background = np.array([])
means_data, meds_data = [], []
for i in G2_arr:
    means_data.append(np.mean(i))
    meds_data.append(np.median(i))
    background = np.concatenate((background, bootstrap_resample(i, n_min)))

n = 10000
matrix = np.array([])
means_boot, meds_boot = [], []
for b in range(n):
    subset = bootstrap_resample(background, n_min)
    matrix = np.concatenate((matrix, subset))
    means_boot.append(np.mean(subset))
    meds_boot.append(np.median(subset))

matrix = matrix.reshape(n, n_min)  # not necessary anymore

means_pooled = np.concatenate((means_boot, means_data))
meds_pooled = np.concatenate((meds_boot, meds_data))

sorted_means = np.sort(means_pooled)
sorted_meds = np.sort(meds_pooled)

p_means, p_medians = [], []
plt.figure(figsize=(12,6))
plt.title("Empirical p-value (calculated with means)")
plt.plot(means_pooled, 'b.', label = "means pooled")
plt.plot(sorted_means, 'r.', label = "means pooled (sorted)")
for n,val in enumerate(means_data):
    x = list(sorted_means).index(val)
    perc = x / len(sorted_means)
    if perc > 0.5:
        perc = (1 - perc) * 2
    else:
        perc *= 2
    plt.plot(x, val, '*', markersize=20, label=str(sub_keys[n])+" (emperical p-value = "+str(round(perc,2))+")")
    p_means.append(perc)
plt.xlabel("Sample number")
plt.ylabel("G2 (kT)")
plt.legend(prop={'size': 10})
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
plt.savefig(folder+"Emperical p-value (means)",dpi=600)
# plt.show()
plt.close()

plt.figure(figsize=(12,6))
plt.title("Empirical p-value (calculated with medians)")
plt.plot(meds_pooled, 'b.', label = "meds pooled")
plt.plot(sorted_meds, 'r.', label = "meds pooled (sorted)")
for n,val in enumerate(meds_data):
    x = list(sorted_meds).index(val)
    perc = x / len(sorted_meds)
    if perc > 0.5:
        perc = (1 - perc) * 2
    else:
        perc *= 2
    plt.plot(x, val, '*', markersize=20, label=str(sub_keys[n])+" (emperical p-value = "+str(round(perc,2))+")")
    p_medians.append(perc)
plt.legend(prop={'size': 10})
plt.savefig(folder+"Emperical p-value (median)",dpi=600)
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
plt.xlabel("Sample number")
plt.ylabel("G2 (kT)")
# plt.show()
plt.close()

plt.figure(figsize=(15,10))
plt.plot(NRLs, p_means, '-o', label="Emperical p-value (mean)")
plt.plot(NRLs, p_medians, '-o',  label="Emperical p-value (median)")
plt.legend(frameon=False)
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
# plt.show()
plt.ylabel("Emperical p-value")
plt.ylim(-0.1,1.1)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.xticks(int_keys, sub_keys, rotation=45)
plt.savefig(folder+"Emperical p-value scatter-plot",dpi=600)
# plt.show()
# plt.close()

