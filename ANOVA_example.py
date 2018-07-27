import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

folder = "C:\\Users\\tbrouwer\\Desktop\\clean\\"
df_G2 = pd.read_csv(folder+"G2_pars_all.dat", sep='\t')

data = []
keys = list(df_G2.keys())
for n,key in enumerate(keys):
    sub_data = df_G2[key]
    sub_data = np.array(sub_data[~np.isnan(sub_data)])
    data.append(sub_data)

boxplot, sub_keys, int_keys = [], [], []
NRLs = [3,4,6,8,10,11,12]
for n,i in enumerate(NRLs):
    # print(keys[i],data[i])
    boxplot.append(data[i])
    sub_keys.append(keys[i])
    NRL = int(keys[i][:3])
    int_keys.append(NRL)

    statistic, norm_p_val = stats.normaltest(data[i])
    print(str(keys[i])+" - normaltest value  =", norm_p_val)

    x = np.random.normal(NRL, 0.05, size=len(data[i]))
    plt.plot(x, data[i], '.', color="black", zorder=10, alpha=0.2)


plt.boxplot(boxplot, showfliers=False, positions=int_keys)

plt.xticks(int_keys,sub_keys, rotation=45)
plt.title("G2 (kT)", fontsize=20)

f_val, p_val = stats.f_oneway(boxplot[0],boxplot[1],boxplot[2],boxplot[3],boxplot[4],boxplot[5],boxplot[6])
print("One-way ANOVA P =", p_val)

plt.show()

