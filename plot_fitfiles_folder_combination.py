import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import numpy as np
import matplotlib.cm as cm

fitfile_path = "N:\\Brouwer\\NRL\\assembly\\Figure 1\\"
# fitfile_path = "N:\\Brouwer\\NRL\\best\\combination\\figure_2B\\"
save_path = fitfile_path

standard_trajectory = True
evaluate_ruptures = False

fitfiles = []
os.chdir(fitfile_path)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

colors = iter(cm.rainbow(np.linspace(0, 1, len(fitfiles))))
plt.figure(figsize=(30, 20))

for n, fitfile in enumerate(fitfiles):

    p = {}

    print("Processing fitfile... " + str(fitfile))

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile, p)

    # open data
    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_plain(fitfile_path, fitfile, p, standard_trajectory=standard_trajectory, evaluate_ruptures=evaluate_ruptures)
    wlc = np.transpose(transitions[2])[-1]

    # offset all curves to L_c 177
    offset = (177 - p['NRL']) * 16 * 0.33/1000

    # offset all curves to z[f=10]
    for n1,z in enumerate(z_fit_pull):
        if f_pull[n1] > 11:
            offset = -(z - 1.2)
            break

    # no offset
    # offset = 0
    offset -= (0.02 * n)

    # correcting the offset
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
    plt.xlim(0.25,1.25)

    if p['NRL'] == 167:
        plt.scatter(z_pull, f_pull, color='red', label=str(int(p['NRL'])), s=30, zorder=25,
                    facecolors='none')
    elif p['NRL'] == 176:
        # plt.scatter(z_pull, f_pull, color='firebrick', label=str(int(p['NRL']))+" "+ fitfile[:-4], s=30, zorder=25,
        #             facecolors='none')
        plt.scatter(z_pull, f_pull, color=next(colors), label=str(int(p['NRL'])) + " " + fitfile[:-4], s=30, zorder=25,
                    facecolors='none')
    elif p['NRL'] == 177:
        # plt.scatter(z_pull, f_pull, color='salmon', label=str(int(p['NRL'])), s=30, zorder=25,
        #             facecolors='none')
        plt.scatter(z_pull, f_pull, color=next(colors), label=str(int(p['NRL'])) +" "+ fitfile[:-4], s=30, zorder=25,
                    facecolors='none')
    elif p['NRL'] == 168:
        plt.scatter(z_pull, f_pull, color='grey', label=str(int(p['NRL'])), s=30, zorder=25,
                    facecolors='none')
    else:
        plt.scatter(z_pull, f_pull, color=next(colors), label=str(int(p['NRL']))+" "+fitfile[:-4], s=30, zorder=25, facecolors='none')
        # plt.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none')

    if fitfile == '2018_08_08 177x16_data_002_7.fit':
        f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles_plain(fitfile_path, fitfile, p,
                                                                                               standard_trajectory=False,
                                                                                               evaluate_ruptures=evaluate_ruptures)
        offset = 00
        z_pull -= offset
        plt.scatter(z_pull, f_pull, color='red', label=str(int(p['NRL'])) + fitfile[:-4], s=30, zorder=25,
                    facecolors='none')
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)
    # plt.plot(wlc, f_pull, '--', color="black", zorder=100)
    plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)

    # select highest number of steps for transitions
    if n == 0:
        # plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        # # plot transitions
        # for t in range(len(np.transpose(transitions[0]))):
        #     # plt.plot(np.transpose(transitions[0])[t]+offset,f_trans,'--',color='lightgrey')  # transition 1
        #     # plt.plot(np.transpose(transitions[1])[t]+offset,f_trans,'--',color='lightgrey')  # transition 2
        #     plt.plot(np.transpose(transitions[2])[t] + offset, f_pull, '--', color='lightgrey')  # transition 3

plt.title("NRL combination")
plt.legend(loc=2, frameon=False)
plt.savefig(save_path + "NRL combination")
plt.show()