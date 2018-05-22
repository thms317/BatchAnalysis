import os
import matplotlib.pyplot as plt
import functions as func
import ba_tools as ba
import glob
import numpy as np


def default_pars():
    pars = {}
    pars['kT'] = 4.114  # pN nm
    pars['L0'] = 0.34  # nm / base pair
    pars['L_bp'] = 4753   # number of base pairs
    pars['P_nm'] = 50  # persistence length
    pars['S_pN'] = 1000  # stretch modulus
    pars['z0_nm'] = 0  # offset in nm / subunit
    pars['NRL'] = "168x16"  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    pars['global_drift'] = True
    return pars

p = default_pars()

def main():
    plt.close("all")

    table_path = "C:\\Users\\brouw\\Desktop\\Data\\"
    table_file = "180509"

    measurements = ba.build_measurements(table_path + table_file + ".xlsx", p)

    for measurement in measurements:
        f_pull, z_pull, f_release, z_release, title = ba.read_analyze(measurement, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])

        plt.ylabel('F (pN)')
        plt.xlabel('z (nm)')
        plt.tick_params(direction='in', top=True, right=True)
        plt.ylim(-1, 65)
        plt.xlim(500, 2300)

        plt.scatter(1000 * z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        plt.scatter(1000 * z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        plt.plot(wlc, f_wlc, '--', color="black", label="WLC", zorder=10000)
        plt.plot([], [], ' ', label="Drift: " + str(p['drift']))  # quick and very dirty

        plt.legend(loc=2, frameon=False)
        plt.title(title)

        if p['save'] == True:
            new_path = table_path + table_file + "\\Figs\\"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            plt.savefig(new_path + title)

        # plt.show()
        plt.close()

    return

def main_rot():
    plt.close("all")

    table_path = "C:\\Users\\brouw\\Desktop\\Data\\"
    table_file = "180509"

    measurements = ba.build_measurements(table_path + table_file + ".xlsx", p)

    for measurement in measurements:

        f_pull, z_pull, f_release, z_release, title = ba.read_analyze(measurement, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])
        drift_FE = p['drift']

        twist_pos, twist_neg, z_pos, z_neg, lnd_pos, lnd_neg = ba.read_analyze_rot(measurement, p)
        drift_rot = p['drift']

        fig = plt.figure(figsize=(15,5))

        # rotation
        ax1 = fig.add_subplot(1,2,1)

        ax1.set_ylabel('z (nm)')
        ax1.set_xlabel('$\sigma$')
        ax1.tick_params(direction='in', top=True, right=True)
        ax1.set_xlim(-0.075, 0.075)
        ax1.set_ylim(-500, 500)

        # twist on x-axis
        # ax1.scatter(twist_pos, 1000*z_pos, color='darkgreen', label="Forward twisting", s=10, zorder=25, facecolors='none')
        # ax1.scatter(twist_neg, 1000*z_neg, color='lightgrey', label='Backward twisting', s=10, zorder=15, facecolors='none')
        # LND on x-axis
        ax1.scatter(lnd_pos, 1000*z_pos, color='darkgreen', label="Forward twisting", s=10, zorder=25, facecolors='none')
        ax1.scatter(lnd_neg, 1000*z_neg, color='lightgrey', label='Backward twisting', s=10, zorder=15, facecolors='none')

        ax1.plot([], [], ' ', label="Drift: " + str(drift_rot))  # quick and very dirty

        ax1.legend(loc=2, frameon=False)
        ax1.set_title("Twist @ 1 pN")

        # force-extension
        ax2 = fig.add_subplot(1, 2, 2)

        ax2.set_ylabel('F (pN)')
        ax2.set_xlabel('z (nm)')
        ax2.tick_params(direction='in', top=True, right=True)
        ax2.set_ylim(-1, 65)
        ax2.set_xlim(500, 2300)

        ax2.plot([], [], ' ', label="Drift: " + str(drift_FE))  # quick and very dirty
        ax2.scatter(1000 * z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        ax2.scatter(1000 * z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        ax2.plot(wlc, f_wlc, '--', color="black", label="WLC", zorder=10000)

        ax2.legend(loc=2, frameon=False)
        ax2.set_title("Force Extension")

        fig.suptitle(title)

        if p['save'] == True:
            new_path = table_path + table_file + "\\Selected Figs\\"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            fig.savefig(new_path + title)

        # fig.show()
        plt.close()

    # TODO fix why pycharm closes figures

    return


def main_fitfiles():

    fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\Fitfiles\\"

    fitfiles = []
    os.chdir(fitfile_path)
    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    for fitfile in fitfiles:

        f_pull, f_release, z_pull, z_release, z_fit_pull, transitions, f_trans = ba.read_fitfiles(fitfile_path, fitfile, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])

        # plot transitions
        for t in range(len(np.transpose(transitions[0]))):
            # plt.plot(np.transpose(transitions[0])[t],f_trans,'--',color='lightgrey')  # transition 1
            # plt.plot(np.transpose(transitions[1])[t],f_trans,'--',color='lightgrey')  # transition 2
            plt.plot(np.transpose(transitions[2])[t],f_trans,'--',color='lightgrey')  # transition 3

        plt.ylabel('F (pN)')
        plt.xlabel('z (nm)')
        plt.tick_params(direction='in', top=True, right=True)

        plt.ylim(-1, 60)
        plt.xlim(0, 1.8)

        plt.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        plt.scatter(z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        plt.plot(wlc/1000, f_wlc, '--', color="black", label="WLC", zorder=100)
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        plt.legend(loc=2, frameon=False)
        plt.title(fitfile[:-4]+" - fitfile (mind offset!)")

        # TODO fix the offset in title

        if p['save'] == True:
            new_path = fitfile_path + "\\Figs\\"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            plt.savefig(new_path + fitfile[:-4])

        # plt.show()
        plt.close()

    return

def main_assemble_pars():

    logfile_path = "C:\\Users\\brouw\\Desktop\\Data\\Fitfiles\\"

    logfiles = []
    os.chdir(logfile_path)
    for file in glob.glob("*.fit"):
        logfiles.append(file)

    for logfile in logfiles:
        print(logfile)

    # TODO read pars from logfile

if __name__ == "__main__":
    # main()
    main_rot()
    # main_fitfiles()
    # main_assemble_pars()
