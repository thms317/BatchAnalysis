import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import numpy as np
import matplotlib.cm as cm

fitfile_path = "N:\\Brouwer\\NRL\\assembly\\Figure 1\\"
save_path = fitfile_path + "figs\\"
if not os.path.exists(save_path):
    os.makedirs(save_path)

standard_trajectory = False
evaluate_ruptures = False
individual_figures = False

fitfiles = []
os.chdir(fitfile_path)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

colors = iter(cm.rainbow(np.linspace(0, 1, len(fitfiles))))

# if assembly
plt.figure(figsize=(30, 20))

for n, fitfile in enumerate(fitfiles):

    if individual_figures:
        # if separate figures
        plt.figure(figsize=(15, 10))

    p = {}

    print("Processing fitfile... " + str(fitfile))

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile, p)

    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_assembly(fitfile_path, fitfile, p, standard_trajectory=standard_trajectory, evaluate_ruptures=evaluate_ruptures)
    wlc = np.transpose(transitions[2])[-1]

    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    plt.ylabel('F (pN)')
    plt.xlabel('z (nm)')
    plt.tick_params(direction='in', top=True, right=True)

    plt.ylim(-1, 60)
    plt.xlim(0,2)

    if individual_figures:
        # plot transitions
        for t in range(len(np.transpose(transitions[0]))):
            # plt.plot(np.transpose(transitions[0])[t],f_trans,'--',color='lightgrey')  # transition 1
            # plt.plot(np.transpose(transitions[1])[t],f_trans,'--',color='lightgrey')  # transition 2
            plt.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3

    plt.scatter(z_pull, f_pull, color=next(colors), label=str(int(p['NRL'])), s=30, zorder=25, facecolors='none')
    # plt.scatter(z_pull, f_pull, color='red', label=str(int(p['NRL'])), s=30, zorder=25, facecolors='none')
    # plt.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none')
    # plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
    if individual_figures:
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)
    else:
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)

    if individual_figures:
            # print pars in figure
        report = str(table[0]) + '\n' + str(table[1]) + '\n' + str(table[6]) + '\n' + str(table[2]) + '\n' + str(table[3]) + '\n' + str(
            table[4]) + '\n' + str(table[5]) + '\n' + str(table[7])
        plt.annotate(report, xy=(0, 0.75), xytext=(12, -12), va='top', xycoords='axes fraction',
                     textcoords='offset points')

        plt.title("Fitfile: " + fitfile)
        plt.legend(loc=2, frameon=False)
        plt.savefig(save_path + str(int(p['NRL'])), dpi=600)
        plt.close()

if individual_figures == False:
    plt.title("NRL combination")
    plt.legend(loc=2, frameon=False)
    plt.savefig(save_path + "combination", dpi=600)
    plt.close()