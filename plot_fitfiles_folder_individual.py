import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import numpy as np
import matplotlib.cm as cm

fitfile_path = "C:\\Users\\brouw\\Desktop\\NRL\\167x16\\"
fitfile_path = "C:\\Users\\brouw\\Desktop\\NRL\\168x16\\"
fitfile_path = "C:\\Users\\brouw\\Desktop\\NRL\\176x16\\"
fitfile_path = "C:\\Users\\brouw\\Desktop\\NRL\\177x16\\"
save_path = fitfile_path + "figs\\"
if not os.path.exists(save_path):
    os.makedirs(save_path)

standard_trajectory = False
evaluate_ruptures = False

fitfiles = []
os.chdir(fitfile_path)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

color = "darkgreen"
color = "navy"
color = "red"
color = "orange"

for n, fitfile in enumerate(fitfiles):

    plt.figure(figsize=(30, 20))

    p = {}

    print("Processing fitfile... " + str(fitfile))

    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_plain(fitfile_path, fitfile, standard_trajectory=standard_trajectory, evaluate_ruptures=evaluate_ruptures)
    wlc = np.transpose(transitions[2])[-1]

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile, p)

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

    plt.scatter(z_pull, f_pull, color=color, label=str(int(p['NRL'])), s=30, zorder=25, facecolors='none')
    # plt.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none')
    plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
    plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

    # print pars in figure
    report = str(table[0]) + '\n' + str(table[1]) + '\n' + str(table[6]) + '\n' + str(table[2]) + '\n' + str(table[3]) + '\n' + str(
        table[4]) + '\n' + str(table[5]) + '\n' + str(table[7])
    plt.annotate(report, xy=(0, 0.75), xytext=(12, -12), va='top', xycoords='axes fraction',
                 textcoords='offset points')

    plt.title(fitfile[:-4])
    plt.legend(loc=2, frameon=False)
    plt.savefig(save_path + fitfile[:-4])
    plt.close()