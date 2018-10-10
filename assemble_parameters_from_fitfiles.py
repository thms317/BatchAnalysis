import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import pandas as pd

p = {}

plt.close("all")

assemble_pars_file = False

save_folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative Histograms\\"

fig = plt.figure(figsize=(30, 18))
plt.rcParams.update({'font.size': 20})
plt.rc('axes', linewidth=3)

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\195x16\\"
p['NRL_str'] = "195x16"

fitfiles = []
os.chdir(folder)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

ass_fit_pars, ass_fit_errors, ass_dict = [], [], []

for fitfile in fitfiles:
    print("Processing fitfile... " + str(fitfile))

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors = ba.read_logfile_clean(folder, logfile, p)

    ass_fit_pars.append(fit_pars)
    ass_fit_errors.append(fit_errors)
    ass_dict.append(p)

# assemble parameters into histogram
ba.plot_combined_hist(fig, ass_fit_pars, ass_fit_errors, folder, p, color='blue',zorder=100, ymax=30)

if assemble_pars_file:
    values = []
    for n, dict in enumerate(ass_dict):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(folder + p['NRL_str']+ "_assembled_pars.dat", sep='\t')


# save histograms
# fig.suptitle("Mechanical Parameters of Chromatin are independent of Numbers of Repeats \n 197x15: n = "+str(n_15)+"    |    197x25: n = "+str(n_25))
# fig.suptitle("Mechanical Parameters of Chromatin are independent of Numbers of Repeats")
plt.suptitle(p['NRL_str'] + " (n = "+ str(len(fitfiles)) + ")")
plt.savefig(save_folder+p['NRL_str']+"_combined_pars")
