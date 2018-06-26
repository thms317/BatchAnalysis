import os
import matplotlib.pyplot as plt
import functions as func
import ba_tools as ba
import glob
import numpy as np
import csv

def default_pars():
    pars = {}
    pars['kT'] = 4.114  # pN nm
    pars['L0'] = 0.34  # nm / base pair
    pars['L_bp'] = 4769  # number of base pairs
    pars['P_nm'] = 50  # persistence length
    pars['S_pN'] = 1000  # stretch modulus
    pars['z0_nm'] = 0  # offset in nm / subunit
    pars['NRL'] = 169  # nucleosome repeat length
    pars['repeats'] = 16  # number of repeats
    pars['type'] = "Human"  # type of histone
    pars['NRL_str'] = str(pars['NRL'])+'x'+str(pars['repeats'])+'_'+pars['type']  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    return pars

p = default_pars()

def main_measurement_files():
    plt.close("all")

    table_path = "C:\\Users\\brouw\\Desktop\\Data\\"
    table_file = "180622_171"

    measurements = ba.build_measurements(table_path, table_file + ".txt", p)
    drift_arr = []
    z_offset_arr = []
    for measurement in measurements:
        print("Processing measurement... " + str(measurement))
        f_pull, z_pull, f_release, z_release, title, drift, z_offset = ba.read_analyze(measurement, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])
        drift_FE = p['drift']
        drift_arr.append(drift)
        z_offset_arr.append(z_offset)

        twist_pos, twist_neg, z_pos, z_neg, lnd_pos, lnd_neg = ba.read_analyze_rot(measurement, p)
        drift_rot = p['drift']

        fig = plt.figure(figsize=(15, 5))

        # rotation
        ax1 = fig.add_subplot(1, 2, 1)

        ax1.set_ylabel('z (nm)')
        ax1.set_xlabel('$\sigma$')
        ax1.tick_params(direction='in', top=True, right=True)
        ax1.set_xlim(-0.025, 0.025)
        ax1.set_ylim(-500, 500)

        # twist on x-axis
        # ax1.scatter(twist_pos, 1000*z_pos, color='darkgreen', label="Forward twisting", s=10, zorder=25, facecolors='none')
        # ax1.scatter(twist_neg, 1000*z_neg, color='lightgrey', label='Backward twisting', s=10, zorder=15, facecolors='none')
        # LND on x-axis
        ax1.scatter(lnd_pos, 1000 * z_pos, color='darkgreen', label="Forward twisting", s=10, zorder=25,
                    facecolors='none')
        ax1.scatter(lnd_neg, 1000 * z_neg, color='lightgrey', label='Backward twisting', s=10, zorder=15,
                    facecolors='none')

        ax1.plot([], [], ' ', label="Drift: " + str(drift_rot))  # quick and very dirty

        ax1.legend(loc=2, frameon=False)
        ax1.set_title("Twist @ 1 pN")

        # force-extension
        ax2 = fig.add_subplot(1, 2, 2)

        ax2.set_ylabel('F (pN)')
        ax2.set_xlabel('z (nm)')
        ax2.tick_params(direction='in', top=True, right=True)
        ax2.set_ylim(-1, 65)
        ax2.set_xlim(0, 2000)

        ax2.plot([], [], ' ', label="Drift: " + str(drift_FE))  # quick and very dirty
        ax2.scatter(1000 * z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        ax2.scatter(1000 * z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        ax2.plot(wlc, f_wlc, '--', color="black", label="WLC", zorder=10000)

        ax2.legend(loc=2, frameon=False)
        ax2.set_title("Force Extension")

        fig.suptitle(title)

        if p['save'] == True:
            new_path = table_path + table_file[:6] + "\\Selected figures" + table_file[6:] + "\\"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            fig.savefig(new_path + title)

        # fig.show()
        plt.close()

    if p['save'] == True:
        measurements = np.transpose(measurements)
        measurements = np.append(measurements, np.full((len(z_offset_arr)), p['L_bp']))  # append length
        measurements = np.append(measurements, np.full((len(z_offset_arr)), p['NRL']))  # append NRL
        measurements = np.append(measurements, np.full((len(z_offset_arr)), p['repeats']))  # append repeats
        measurements = np.append(np.transpose(measurements), z_offset_arr)  # append z-offset
        measurements = np.transpose(np.append(measurements, drift_arr).reshape((9, -1)))  # append drift
        headers = ["date", "data_", "bead", "type", "length (bp)", "NRL", "repeats", "z-offset", "drift"]

        measurements = np.append(headers,measurements).reshape((-1, 9))
        with open(new_path + table_file + "_list_of_measurements.txt", "w+") as my_csv:
            csvWriter = csv.writer(my_csv, delimiter='\t')
            csvWriter.writerows(measurements)

    return


def main_fitfiles():
    fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\180621\\Fitfiles (refitting G2)\\"

    fitfiles = []
    os.chdir(fitfile_path)
    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    ass_fit_pars = []
    ass_fit_errors = []

    for fitfile in fitfiles:
        print("Processing fitfile... " + str(fitfile))

        f_pull, f_release, z_pull, z_release, z_fit_pull, transitions = ba.read_fitfiles(fitfile_path, fitfile, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(55.2896842957)), 1000)
        # f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)  # this one is better, although some measurements crash
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])

        # read pars from logfile
        logfile = fitfile[:-3] + "log"
        fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile)
        ass_fit_pars.append(fit_pars)
        ass_fit_errors.append(fit_errors)

        fig = plt.figure(figsize=(30, 10))
        plt.rcParams.update({'font.size': 20})
        plt.rc('axes', linewidth=3)

        # number of nucleosomes
        ax0 = fig.add_subplot(1, 2, 1)

        ax0.set_ylabel('F (pN)')
        ax0.set_xlabel('z (nm)')
        ax0.tick_params(direction='in', top=True, right=True)

        ax0.set_ylim(0, 6)
        ax0.set_xlim(0.25, 1.25)

        ax0.set_title("Zoom in")

        ax0.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=30, zorder=25, facecolors='none')
        ax0.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, label='Release', facecolors='none')
        ax0.plot(wlc / 1000, f_wlc, '--', color="black", label="WLC", zorder=100)
        ax0.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        ax0.legend(loc=2, frameon=False)

        # number of tetrasomes
        ax1 = fig.add_subplot(1, 2, 2)

        ax1.set_ylabel('F (pN)')
        ax1.set_xlabel('z (nm)')
        ax1.tick_params(direction='in', top=True, right=True)

        ax1.set_ylim(-1, 60)
        ax1.set_xlim(0, 1.8)

        ax1.set_title("Zoom out")

        # print pars in figure
        report = str(table[0]) + '\n' + str(table[1]) + '\n' + str(table[2]) + '\n' + str(table[3]) + '\n' + str(
            table[4]) + '\n' + str(table[5])
        ax1.annotate(report, xy=(0, 0.75), xytext=(12, -12), va='top', xycoords='axes fraction',
                     textcoords='offset points')

        # plot transitions
        for t in range(len(np.transpose(transitions[0]))):
            # ax1.plot(np.transpose(transitions[0])[t],f_trans,'--',color='lightgrey')  # transition 1
            # ax1.plot(np.transpose(transitions[1])[t],f_trans,'--',color='lightgrey')  # transition 2
            ax1.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3

        ax1.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=30, zorder=25, facecolors='none')
        ax1.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, label='Release', facecolors='none')
        ax1.plot(wlc / 1000, f_wlc, '--', color="black", label="WLC", zorder=100)
        ax1.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        ax1.legend(loc=2, frameon=False)

        fig.suptitle(fitfile[:-4] + " - " + p['NRL_str'])

        if p['save'] == True:
            new_path = fitfile_path + "\\Fitfile figures\\"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            plt.savefig(new_path + fitfile[:-4])

        # plt.show()
        plt.close()

    # assemble parameters into histogram
    if p['save'] == True:
        ba.plot_hist(ass_fit_pars, ass_fit_errors, new_path, p, show_plot=False)

    return


if __name__ == "__main__":
    # main_measurement_files()
    main_fitfiles()

