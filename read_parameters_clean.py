import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters\\clean\\"
save_folder = "C:\\Users\\brouw\\Desktop\\Boxplots\\"

df_k = pd.read_csv(folder+"k_pars_all.dat", sep='\t')
df_G1 = pd.read_csv(folder+"G1_pars_all.dat", sep='\t')
df_G2 = pd.read_csv(folder+"G2_pars_all.dat", sep='\t')

keys = list(df_k.keys())

# stiffness

boxplot = []
for n,key in enumerate(keys):
    data = df_k[key]
    data = np.array(data[~np.isnan(data)])
    boxplot.append(data)

    # histograms

    plt.ylabel('count',fontsize=20)
    plt.xlabel('Fiber stiffness (pN/nm)',fontsize=20)
    plt.tick_params(direction='in', top=True, right=True)
    plt.title(key+" (n = "+str(len(data))+")", fontsize=20)
    binwidth = 0.1  # pN/nm
    plt.xlim(0, 2)
    plt.ylim(0, 25)
    try:
        plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth),
             edgecolor='black', color='red', linewidth=1.2)
        plt.savefig(save_folder+key+"_stiffness.png")
        plt.close()
    except:
        plt.close()

plt.figure(figsize=(12,6))
plt.boxplot(boxplot)
plt.ylim(0.1,3)
plt.semilogy()
plt.xticks(np.arange(len(keys))+1,keys, rotation=45)
plt.title("Stiffness (pN/nm)", fontsize=20)

plt.savefig(save_folder+"k_boxplot.png")

# plt.show()
plt.close()

# G1

boxplot = []
for n,key in enumerate(keys):
    data = df_G1[key]
    data = np.array(data[~np.isnan(data)])
    boxplot.append(data)

    # histograms

    binwidth = 1  # kT
    plt.ylabel('count',fontsize=20)
    plt.xlabel('G1 (kT)',fontsize=20)
    plt.tick_params(direction='in', top=True, right=True)
    plt.title(key+" (n = "+str(len(data))+")", fontsize=20)

    plt.xlim(0, 30)
    plt.ylim(0, 25)
    try:
        plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth),
             edgecolor='black', color='green', linewidth=1.2)
        plt.savefig(save_folder+key+"_G1.png")
        plt.close()
    except:
        plt.close()

plt.figure(figsize=(12,6))
plt.boxplot(boxplot)
# plt.ylim(0.1,3)
# plt.semilogy()
plt.xticks(np.arange(len(keys))+1,keys, rotation=45)
plt.title("G1 (kT)",fontsize=20)

plt.savefig(save_folder+"G1_boxplot.png")

# plt.show()
plt.close()

# G2

boxplot = []
for n,key in enumerate(keys):
    data = df_G2[key]
    data = np.array(data[~np.isnan(data)])
    boxplot.append(data)

    # histograms

    binwidth = 1  # kT
    plt.ylabel('count',fontsize=20)
    plt.xlabel('G2 (kT)',fontsize=20)
    plt.tick_params(direction='in', top=True, right=True)
    plt.title(key+" (n = "+str(len(data))+")", fontsize=20)

    plt.xlim(0, 30)
    plt.ylim(0, 25)
    try:
        plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth),
             edgecolor='black', color='blue', linewidth=1.2)
        plt.savefig(save_folder+key+"_G2.png")
        plt.close()
    except:
        plt.close()

plt.figure(figsize=(12,6))
plt.boxplot(boxplot)
# plt.ylim(0.1,3)
# plt.semilogy()
plt.xticks(np.arange(len(keys))+1,keys, rotation=45)
plt.title("G2 (kT)", fontsize=20)

plt.savefig(save_folder+"G2_boxplot.png")

# plt.show()
plt.close()