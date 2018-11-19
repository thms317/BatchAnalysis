import matplotlib.pyplot as plt
import ba_tools as ba
import numpy as np
import measurements_WT_tailless_197 as meas

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\"
save_path = "C:\\Users\\brouw\\Desktop\\WT tailless comparison\\"

evaluate_ruptures = False

measurements = meas.measurements()

# colors = iter(cm.rainbow(np.linspace(0, 1, len(measurements))))
plt.figure(figsize=(15, 10))

for n, measurement in enumerate(measurements):

    fitfile = measurement[3] + "\\" + measurement[0] + "_data_" + measurement[1] + "_" + measurement[2] + ".fit"
    fitfile_ID = measurement[3] + "_" + measurement[0] + "_data_" + measurement[1] + "_" + measurement[2]

    p = {}

    print("Processing fitfile... " + str(measurement))

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(folder, logfile, p)

    # open data
    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_plain(folder, fitfile, p, standard_trajectory=True, evaluate_ruptures=evaluate_ruptures)
    wlc = np.transpose(transitions[2])[-1]

    offset = 0

    # offset all curves to z[f=10]
    for n1, z in enumerate(z_fit_pull):
        if f_pull[n1] > 10:
            offset = -(z - 1.2)
            break

    # # correcting the offset
    z_pull += offset
    z_release += offset
    z_fit_pull += offset
    wlc += offset

    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    plt.ylabel('F (pN)')
    plt.xlabel('z (nm)')
    plt.tick_params(direction='in', top=True, right=True)

    plt.ylim(0, 10)
    plt.xlim(0.55,1.25)

    plt.scatter(z_pull, f_pull, label=fitfile_ID, s=30, zorder=25)
    # if n == 0:
    #     plt.scatter(z_pull, f_pull, label=fitfile_ID, s=30, zorder=25, color='navy', facecolors='none')
    # if n == 1:
    #     plt.scatter(z_pull, f_pull, label=fitfile_ID, s=30, zorder=25, color='red', facecolors='none')

    plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)

    # plotting for the legend
    if n == 0:
        # plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        # # plot transitions - select highest number of steps for transitions
        # for t in range(len(np.transpose(transitions[0]))):
        #     # plt.plot(np.transpose(transitions[0])[t]+offset,f_trans,'--',color='lightgrey')  # transition 1
        #     # plt.plot(np.transpose(transitions[1])[t]+offset,f_trans,'--',color='lightgrey')  # transition 2
        #     plt.plot(np.transpose(transitions[2])[t] + offset, f_pull, '--', color='lightgrey')  # transition 3

plt.title(measurement[3][:6] + " WT vs. tailless")
plt.legend(loc=2, frameon=False)
plt.savefig(save_path + measurement[3][:6]+ " WT vs tailless (combination)", dpi=600)
plt.show()