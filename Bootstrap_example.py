import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters\\clean\\"
save_folder = "C:\\Users\\tbrouwer\\Desktop\\Boxplots\\"

df_G2 = pd.read_csv(folder + "G2_pars_all.dat", sep='\t')

data = []
keys = list(df_G2.keys())
for n, key in enumerate(keys):
    sub_data = df_G2[key]
    sub_data = np.array(sub_data[~np.isnan(sub_data)])
    data.append(sub_data)

boxplot, sub_keys, int_keys, length = [], [], [], []
G2_tot = np.array([])
NRLs = [3, 4, 6, 8, 10, 12]
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
    # plt.plot(x, data[i], '.', color="black", zorder=10, alpha=0.2)

# plt.boxplot(boxplot, showfliers=False, positions=int_keys)
# plt.xticks(int_keys,sub_keys, rotation=45)
# plt.title("G2 (kT)", fontsize=20)
# plt.show()

# remove lowest number of measurements
# length.remove(min(length))
n_min = min(length)
print(length)
def bootstrap_resample(X, n):
    X = np.array(X)
    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    X_resample = np.array(X[resample_i])
    return X_resample

background = np.array([])
means_data, meds_data = [], []
for i in boxplot:
    means_data.append(np.mean(i))
    meds_data.append(np.median(i))
    background = np.concatenate((background, bootstrap_resample(i, n_min)))

n = 1000
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
plt.xlabel("Sample number")
plt.ylabel("G2 (kT)")
plt.legend()
plt.savefig(save_folder+"Emperical p-value (means)",dpi=600)
plt.show()
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
plt.legend()
plt.savefig(save_folder+"Emperical p-value (median)",dpi=600)
plt.xlabel("Sample number")
plt.ylabel("G2 (kT)")
# plt.show()
plt.close()