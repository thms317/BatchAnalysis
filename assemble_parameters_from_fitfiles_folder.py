import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import pandas as pd
import numpy as np

assemble_pars_file = True
analysis = False

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\"
save_folder = "C:\\Users\\tbrouwer\\Desktop\\Chromatin Force Spectroscopy\\"

# fig = plt.figure(figsize=(30, 18))
# plt.rcParams.update({'font.size': 20})
# plt.rc('axes', linewidth=3)

series = np.arange(167,178)

global_fit_pars, global_fit_errors, ass_dict = [], [], []

for NRL in series:
    print("Processing NRL... " + str(NRL))

    subfolder = folder + str(NRL) + "x16\\"

    fitfiles, ass_fit_pars, ass_fit_errors = [], [], []

    try:
        os.chdir(subfolder)
    except:
        pass

    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    for fitfile in fitfiles:
        print("Processing fitfile... " + str(fitfile))
        p = {}
        p['NRL_str'] = str(NRL) + "x16"

        # read pars from logfile
        logfile = fitfile[:-3] + "log"
        fit_pars, fit_errors = ba.read_logfile_clean(subfolder, logfile, p)

        ass_fit_pars.append(fit_pars)
        ass_fit_errors.append(fit_errors)
        ass_dict.append(p)

    global_fit_pars.append(ass_fit_pars)
    global_fit_errors.append(ass_fit_errors)

    # assemble parameters into histogram
    ba.plot_hist(ass_fit_pars, ass_fit_errors, save_folder, p, show_plot = False)


if assemble_pars_file:
    values = []
    for n, dict in enumerate(ass_dict):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(save_folder + "assembled_pars.dat", sep='\t')
