import os, glob
import matplotlib.pyplot as plt
import numpy as np
import ba_tools as ba

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\"
# savefolder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative Boxplots\\"
savefolder = "C:\\Users\\brouw\\Desktop\\NRL analysis\\"

save = True
depict_n = False
reference = True

# actual NRLs
NRLs = ["167x16_WT", "167x16_tailless", "197x16_WT", "197x16_tai    lless"]
NRLs_loc = [0, 1, 2.25, 3.25]

plt.rcParams.update({'font.size': 20})  # legend + title size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

boxplot_k, boxplot_G1, boxplot_G2, int_keys, sub_keys = [], [], [], [], []
stiffness_167, G1_167, G2_167, stiffness_197, G1_197, G2_197 = [], [], [], [], [], []

total = 0

if reference:

    for n, ref in enumerate(["167x16", "197x16"]):

        print("Processing NRL... " + str(ref))

        subfolder = folder + str(ref) + "\\"

        logfiles = []
        os.chdir(subfolder)
        for file in glob.glob("*.log"):
            logfiles.append(file)

        p = {}

        for logfile in logfiles:

            try:
                fit_pars, fit_errors = ba.read_logfile_clean(subfolder, logfile, p)

                if n == 0:
                    stiffness_167.append(fit_pars[2])
                    G1_167.append(fit_pars[3])
                    G2_167.append(fit_pars[4])
                elif n == 1:
                    stiffness_197.append(fit_pars[2])
                    G1_197.append(fit_pars[3])
                    G2_197.append(fit_pars[4])

            except:
                pass

for n, NRL in enumerate(NRLs):

    print("Processing NRL... " + str(NRL))

    subfolder = folder + str(NRL) + "\\"

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

            total += 1

        except:
            pass

    boxplot_k.append(stiffness)
    boxplot_G1.append(G1)
    boxplot_G2.append(G2)

    int_keys.append((NRLs_loc[n]))
    sub_keys.append(str((NRL)))

    x = np.random.normal(NRLs_loc[n], 0.05, size=len(stiffness))
    y = np.random.normal(NRLs_loc[n], 0.05, size=len(G1))
    z = np.random.normal(NRLs_loc[n], 0.05, size=len(G2))

    plt.figure(0, figsize=(12, 7))
    plt.plot(x, stiffness, '.', color="black", zorder=10, alpha=0.2)
    plt.figure(1, figsize=(12, 7))
    plt.plot(y, G1, '.', color="black", zorder=10, alpha=0.2)
    plt.figure(2, figsize=(12, 7))
    plt.plot(z, G2, '.', color="black", zorder=10, alpha=0.2)

# stiffness
plt.figure(0)
plt.boxplot(boxplot_k, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Stiffness (pN/nm)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("Stiffness (kN/nm)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
ymin, ymax = plt.ylim()
plt.vlines(1.625, ymin, ymax, linewidth=3)
if reference:
    xmin, xmax = plt.xlim()
    plt.hlines(np.median(stiffness_167), xmin, 1.625, linestyles='--', label="167x16 human", colors="green")
    plt.hlines(np.median(stiffness_197), 1.625, xmax, linestyles='--', label="197x16 human", colors="red")
    plt.legend(frameon=False, loc='upper center', ncol=2)
if depict_n:
    for n, loc in enumerate(NRLs_loc):
        plt.text(loc, ymax, " n = " + str(len(boxplot_k[n])), size=15, rotation=45, verticalalignment='bottom')
if save:
    plt.savefig(savefolder + "WT_tailless_stiffness_boxplot", dpi=600)
plt.semilogy()
plt.ylim(0.1, 10)
if save:
    plt.savefig(savefolder + "WT_tailless_stiffness_boxplot_log", dpi=600)

# G1
plt.figure(1)
plt.boxplot(boxplot_G1, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Stacking Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("$\Delta G_{1}$ (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
ymin, ymax = plt.ylim()
plt.vlines(1.625, ymin, ymax, linewidth=3)
if reference:
    xmin, xmax = plt.xlim()
    plt.hlines(np.median(G1_167), xmin, 1.625, linestyles='--', label="167x16 human", colors="green")
    plt.hlines(np.median(G1_197), 1.625, xmax, linestyles='--', label="197x16 human", colors="red")
    plt.legend(frameon=False, loc='upper center', ncol=2)
if depict_n:
    for n, loc in enumerate(NRLs_loc):
        plt.text(loc, ymax, " n = " + str(len(boxplot_k[n])), size=15, rotation=45, verticalalignment='bottom')
if save:
    plt.savefig(savefolder + "WT_tailless_G1_boxplot", dpi=600)

# G2
plt.figure(2)
plt.boxplot(boxplot_G2, showfliers=False, positions=int_keys)
plt.xticks(int_keys, sub_keys, rotation=0)
# plt.title("Partial Unwrapping Energy (kT)", fontsize=20)
plt.xlabel("Nucleosome Repeat Length (bp)")
plt.ylabel("$\Delta G_{2}$ (kT)")
plt.tick_params(direction='in', top=True, right=True, length=6, width=3)
ymin, ymax = plt.ylim()
plt.vlines(1.625, ymin, ymax, linewidth=3)
if reference:
    xmin, xmax = plt.xlim()
    plt.hlines(np.median(G2_167), xmin, 1.625, linestyles='--', label="167x16 human", colors="green")
    plt.hlines(np.median(G2_197), 1.625, xmax, linestyles='--', label="197x16 human", colors="red")
    plt.legend(frameon=False, loc='upper center', ncol=2)
if depict_n:
    for n, loc in enumerate(NRLs_loc):
        plt.text(loc, ymax, " n = " + str(len(boxplot_k[n])), size=15, rotation=45, verticalalignment='bottom')
if save:
    plt.savefig(savefolder + "  WT_tailless_G2_boxplot", dpi=600)

print("Total number of measurements... " + str(total))

plt.show()
plt.close()
