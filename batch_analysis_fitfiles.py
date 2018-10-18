import os
import matplotlib.pyplot as plt
import ba_tools as ba
import glob
import numpy as np
import pandas as pd

def default_pars():
    pars = {}
    pars['kT'] = 4.114  # pN nm
    pars['L0'] = 0.34  # nm / base pair
    pars['L_bp'] = 5201  # number of base pairs
    pars['P_nm'] = 50  # persistence length
    pars['S_pN'] = 1000  # stretch modulus
    pars['z0_nm'] = 0  # offset in nm / subunit
    pars['NRL'] = 197  # nucleosome repeat length
    pars['repeats'] = 16  # number of repeats
    pars['type'] = "human"  # type of histone
    pars['NRL_str'] = str(pars['NRL'])+'x'+str(pars['repeats'])+'_'+pars['type']  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    pars['standard'] = True
    pars['radius_um'] = []  # radius of circle (um)
    return pars

p = default_pars()

fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\"

plot_rot = True

M = []

fitfiles = []
os.chdir(fitfile_path)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

ass_fit_pars = []
ass_fit_errors = []

for fitfile in fitfiles:
    m = ba.measurement_pars()

    print("Processing fitfile... " + str(fitfile))

    f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles(fitfile_path, fitfile, p, m)
    wlc = np.transpose(transitions[2])[-1]

    # read pars from logfile
    logfile = fitfile[:-3] + "log"
    fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile, m)
    ass_fit_pars.append(fit_pars)
    ass_fit_errors.append(fit_errors)

    if plot_rot:
        fig = plt.figure(figsize=(30, 20))
    else:
        fig = plt.figure(figsize=(30, 10))

    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    # number of nucleosomes
    if plot_rot:
        ax0 = fig.add_subplot(2, 2, 3)
    else:
        ax0 = fig.add_subplot(1, 2, 1)

    ax0.set_ylabel('F (pN)')
    ax0.set_xlabel('z (nm)')
    ax0.tick_params(direction='in', top=True, right=True)

    ax0.set_ylim(0, 12)
    ax0.set_xlim(0.20, 1.50)

    ax0.set_title("Zoom in")

    ax0.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=30, zorder=25, facecolors='none')
    ax0.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, label='Release', facecolors='none')
    ax0.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

    ax0.legend(loc=2, frameon=False)

    # number of tetrasomes
    if plot_rot:
        ax1 = fig.add_subplot(1, 2, 2)
    else:
        ax1 = fig.add_subplot(1, 2, 2)

    ax1.set_ylabel('F (pN)')
    ax1.set_xlabel('z (nm)')
    ax1.tick_params(direction='in', top=True, right=True)

    ax1.set_ylim(-1, 60)
    ax1.set_xlim(0,2)

    ax1.set_title("Zoom out")

    # print pars in figure
    report = str(table[0]) + '\n' + str(table[1]) + '\n' + str(table[6]) + '\n' + str(table[2]) + '\n' + str(table[3]) + '\n' + str(
        table[4]) + '\n' + str(table[5]) + '\n' + str(table[7])
    ax1.annotate(report, xy=(0, 0.75), xytext=(12, -12), va='top', xycoords='axes fraction',
                 textcoords='offset points')

    # plot transitions
    for t in range(len(np.transpose(transitions[0]))):
        # ax1.plot(np.transpose(transitions[0])[t],f_trans,'--',color='lightgrey')  # transition 1
        # ax1.plot(np.transpose(transitions[1])[t],f_trans,'--',color='lightgrey')  # transition 2
        ax1.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3

    ax1.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=30, zorder=25, facecolors='none')
    ax1.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, label='Release', facecolors='none')
    ax1.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=100)
    ax1.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

    ax1.legend(loc=2, frameon=False)

    if plot_rot:
        measurement = [fitfile[:6],fitfile[12:15],fitfile[16:-4],p['NRL_str']]
        twist_pos, twist_neg, z_pos, z_neg, lnd_pos, lnd_neg = ba.read_analyze_rot(measurement, p)
        drift_rot = p['drift']
        p['radius_um'] = m['radius_um']

        # rotation
        ax2 = fig.add_subplot(2, 2, 1)

        ax2.set_ylabel('z (nm)')
        ax2.set_xlabel('$\sigma$')
        ax2.tick_params(direction='in', top=True, right=True)
        ax2.set_xlim(-0.025, 0.025)
        ax2.set_ylim(-500, 500)

        # twist on x-axis
        ax2.scatter(lnd_pos, 1000 * z_pos, color='darkgreen', label="Forward twisting", s=10, zorder=25,
                    facecolors='none')
        ax2.scatter(lnd_neg, 1000 * z_neg, color='lightgrey', label='Backward twisting', s=10, zorder=15,
                    facecolors='none')

        ax2.plot([], [], ' ', label="Drift: " + str(drift_rot))  # quick and very dirty

        ax2.legend(loc=2, frameon=False)
        ax2.set_title("Twist @ 1 pN")

        fig.suptitle(fitfile[:-4] + " - " + p['NRL_str'])

    if p['save'] == True:
        new_path = fitfile_path + "\\Fitfile figures\\"
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        plt.savefig(new_path + fitfile[:-4])

    # plt.show()
    plt.close()

    M.append(m)

# assemble parameters into histogram
if p['save'] == True:
    ba.plot_hist(ass_fit_pars, ass_fit_errors, new_path, p, show_plot=False)

    values = []
    for n, dict in enumerate(M):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(fitfile_path + p['NRL_str']+ "_assembled_pars.dat", sep='\t')