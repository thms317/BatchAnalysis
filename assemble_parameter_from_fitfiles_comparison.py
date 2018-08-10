import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import pandas as pd
import numpy as np

p = {}

plt.close("all")

assemble_pars_file = True

save_folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters vs Number of Repeats\\"

fig = plt.figure(figsize=(30, 18))  # 4 graphs
# fig = plt.figure(figsize=(15, 9))  # 1 graph
plt.rcParams.update({'font.size': 35})  # legend + title size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=35)
plt.rc('ytick', labelsize=35)


# folder 1

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters vs Number of Repeats\\167 x 16\\"
p['NRL_str'] = "167x16"
# folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters vs Number of Repeats\\197 x 15\\"
# p['NRL_str'] = "197x15"

fitfiles = []
os.chdir(folder)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

n_15 = len(fitfiles)

ass_fit_pars, ass_fit_errors, ass_dict = [], [], []

for fitfile in fitfiles:

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors = ba.read_logfile_clean(folder, logfile, p)

    ass_fit_pars.append(fit_pars)
    ass_fit_errors.append(fit_errors)
    ass_dict.append(p)

# assemble parameters into histogram
ba.plot_combined_hist(fig, ass_fit_pars, ass_fit_errors, folder, p, color='blue',zorder=100)

if assemble_pars_file:
    values = []
    for n, dict in enumerate(ass_dict):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(folder + p['NRL_str']+ "_assembled_pars.dat", sep='\t')


# folder 2

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Parameters vs Number of Repeats\\167 x 30\\"
p['NRL_str'] = "167x30"
# folder = "C:\\Users\\tbrouwer\\Desktop\\Parameters vs Number of Repeats\\197 x 25\\"
# p['NRL_str'] = "197x25"

fitfiles = []
os.chdir(folder)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

n_25 = len(fitfiles)

ass_fit_pars, ass_fit_errors, ass_dict = [], [], []

for fitfile in fitfiles:

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors = ba.read_logfile_clean(folder, logfile, p)

    ass_fit_pars.append(fit_pars)
    ass_fit_errors.append(fit_errors)
    ass_dict.append(p)

# assemble parameters into histogram
ba.plot_combined_hist(fig, ass_fit_pars, ass_fit_errors, folder, p, color='red')

if assemble_pars_file:
    values = []
    for n, dict in enumerate(ass_dict):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(folder + p['NRL_str'] + "_assembled_pars.dat", sep='\t')

# save histograms
# fig.suptitle("Mechanical Parameters of Chromatin are independent of Numbers of Repeats \n 197x15: n = "+str(n_15)+"    |    197x25: n = "+str(n_25))
# fig.suptitle("Mechanical Parameters of Chromatin are independent of Numbers of Repeats")
plt.tight_layout()

# plt.savefig(save_folder+"parameters_197x15_vs_197x25", dpi=600)
plt.savefig(save_folder+"parameters_167x16_vs_167x30", dpi=600)
