import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import numpy as np
import matplotlib.cm as cm

fitfile_path = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\194x16\\"
fitfile_path = "C:\\Users\\brouw\\Desktop\\"
save_path = "C:\\Users\\brouw\\Desktop\\"

standard_trajectory = True
evaluate_ruptures = False

fitfiles = []
os.chdir(fitfile_path)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

colors = iter(cm.rainbow(np.linspace(0, 1, len(fitfiles))))
plt.figure(figsize=(15, 10))

for n, fitfile in enumerate(fitfiles):

    p = {}

    print("Processing fitfile... " + str(fitfile))

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile, p)

    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_plain(fitfile_path, fitfile, p, standard_trajectory=standard_trajectory, evaluate_ruptures=evaluate_ruptures)
    wlc = np.transpose(transitions[2])[-1]

    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    plt.ylabel('F (pN)')
    plt.xlabel('z (nm)')
    plt.tick_params(direction='in', top=True, right=True)

    plt.ylim(-1, 60)
    plt.xlim(0,2)

    # plot transitions
    for t in range(len(np.transpose(transitions[0]))):
        # plt.plot(np.transpose(transitions[0])[t],f_trans,'--',color='lightgrey')  # transition 1
        # plt.plot(np.transpose(transitions[1])[t],f_trans,'--',color='lightgrey')  # transition 2
        plt.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3

    plt.scatter(z_pull, f_pull, color=next(colors), label=fitfile[:-4], s=30, zorder=25, facecolors='none')
    # plt.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none')
    plt.plot(wlc, f_pull, '--', color="black", zorder=100)
    plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)

    if n == 0:
        plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

plt.title("NRL = "+str(int(p['NRL']))+" (assembly)")
plt.legend(loc=2, frameon=False)
plt.savefig(save_path + str(int(p['NRL']))+ " assembly")
plt.close()