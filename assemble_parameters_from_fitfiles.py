import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import pandas as pd
import numpy as np

save_dat = True
save_N = True

# initialize dictionaries
P = []

# initialize figure
fig = plt.figure(figsize=(30, 18))
plt.rcParams.update({'font.size': 20})
plt.rc('axes', linewidth=3)

# overhead
fitfile_folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\167x16_tailless\\"
save_folder = "C:\\Users\\brouw\\Desktop\\"
NRL_string = "167x16_tailless"


fitfiles = []
os.chdir(fitfile_folder)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

ass_fit_pars, ass_fit_errors, ass_dict = [], [], []

for fitfile in fitfiles:
    print("Processing fitfile... " + str(fitfile))

    p = {}
    p['NRL_str'] = NRL_string

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors = ba.read_logfile_clean(fitfile_folder, logfile, p)

    P.append(p)

    ass_fit_pars.append(fit_pars)
    ass_fit_errors.append(fit_errors)
    ass_dict.append(p)

# assemble parameters into histogram
ba.plot_combined_hist(fig, ass_fit_pars, ass_fit_errors, fitfile_folder, p, color='blue',zorder=100, ymax=30)

# save histograms
plt.suptitle(p['NRL_str'] + " (n = "+ str(len(fitfiles)) + ")")
plt.savefig(save_folder+p['NRL_str']+"_histogram")
plt.close()

# save parameters in dat-file
if save_dat:

    values = []
    for n, dict in enumerate(ass_dict):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(save_folder + p['NRL_str']+ "_assembled_pars.dat", sep='\t')

if save_N:

    # initialize figure
    fig = plt.figure(figsize=(30, 9))
    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    ass_fit_pars = np.transpose(np.array(ass_fit_pars))

    # nucleosomes
    ax1 = fig.add_subplot(1, 2, 1)

    ax1.set_ylabel('count')
    ax1.set_xlabel('Nucleosomes in Fiber')
    ax1.tick_params(direction='in', length=6, width=3, top=True, right=True)
    binwidth = 2  # nucleosomes
    ax1.hist(ass_fit_pars[0], bins=np.arange(0, 35, binwidth), edgecolor='black', linewidth=1.2,
             color='red', label=p['NRL_str'], alpha=0.5)

    ax1.legend(loc=1, frameon=False)
    ax1.set_xlim(0, 40)
    ax1.set_ylim(-0.1, 30)
    ax1.yaxis.set_ticks(np.arange(0, 35, 5))

    # tetrasomes
    ax2 = fig.add_subplot(1, 2, 2)

    ax2.set_ylabel('count')
    ax2.set_xlabel('Tetrasomes in Fiber')
    ax2.tick_params(direction='in', length=6, width=3, top=True, right=True)
    binwidth = 2  # tetrasomes
    ax2.hist(ass_fit_pars[1], bins=np.arange(0, 35, binwidth),
             edgecolor='black', linewidth=1.2, color='red', label=p['NRL_str'],alpha=0.5)

    ax2.legend(loc=1, frameon=False)
    ax2.set_xlim(0, 40)
    ax2.set_ylim(-0.1, 30)
    ax2.yaxis.set_ticks(np.arange(0, 35, 5))

    plt.suptitle(p['NRL_str'] + " (n = " + str(len(fitfiles)) + ") \n Nucleosomes vs. Tetrasomes")
    plt.savefig(save_folder + p['NRL_str'] + "_histogram_nuc_tet")
    plt.close()


    plt.scatter(ass_fit_pars[2],ass_fit_pars[3])
    plt.title(NRL_string)
    plt.xlabel("Stiffness (pN/nm)")
    plt.ylabel("$\Delta G{1}$ (kT)")
    plt.tight_layout()
    plt.savefig(save_folder + p['NRL_str'] + "_G1_stiffness")
    plt.close()