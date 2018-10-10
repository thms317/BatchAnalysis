import os, glob
import matplotlib.pyplot as plt
import numpy as np
import ba_tools as ba

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\"
savefolder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative Boxplots\\"
save = True

# actual NRLs
NRLs = list(range(167, 177+1))
NRLs.extend(range(193, 195+1))
# NRL locations (for plotting)
NRLs_loc = list(range(167, 177+1))
NRLs_loc.extend(range(179, 182+1))

plt.rcParams.update({'font.size': 20})  # legend + title size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

boxplot_k, boxplot_G1, boxplot_G2, int_keys, sub_keys = [], [], [], [], []

for n, NRL in enumerate(NRLs):

    print("Processing NRL... " + str(NRL))

    subfolder = folder + str(NRL)+"x16\\"

    logfiles = []
    os.chdir(subfolder)
    for file in glob.glob("*.log"):
        logfiles.append(file)

    stiffness, G1, G2 = [], [], []

    p = {}

    for logfile in logfiles:

        try:
            fit_pars, fit_errors = ba.read_logfile_clean(subfolder, logfile, p)

            stiffness.append(fit_pars[2])
            G1.append(fit_pars[3])
            G2.append(fit_pars[4])
        except:
            pass

    boxplot_k.append(stiffness)
    boxplot_G1.append(G1)
    boxplot_G2.append(G2)

    int_keys.append((int(NRLs_loc[n])))
    sub_keys.append(str(int(NRL)))

    x = np.random.normal(NRLs_loc[n], 0.05, size=len(stiffness))
    y = np.random.normal(NRLs_loc[n], 0.05, size=len(G1))
    z = np.random.normal(NRLs_loc[n], 0.05, size=len(G2))

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
if save:
    plt.savefig(savefolder+"stiffness_boxplot",dpi=600)
plt.semilogy()
plt.ylim(0.1,10)
if save:
    plt.savefig(savefolder+"stiffness_boxplot_log",dpi=600)

plt.figure(1)
plt.boxplot(boxplot_G1, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Stacking Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Stacking Energy (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
if save:
    plt.savefig(savefolder+"G1_boxplot",dpi=600)

plt.figure(2)
plt.boxplot(boxplot_G2, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Partial Unwrapping Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Partial Unwrapping Energy (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
if save:
    plt.savefig(savefolder+"G2_boxplot",dpi=600)

plt.close()
