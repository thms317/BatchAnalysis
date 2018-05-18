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
    pars['NRL'] = "168x16"  # Nucleosome Repeat Length
    pars['drift'] = []
    pars['save'] = True
    return pars


p = default_pars()


def main():
    plt.close("all")

    table_path = "C:\\Users\\brouw\\Desktop\\Data\\"
    table_file = "180518"

    measurements = ba.build_measurements(table_path + table_file + ".xlsx", p)

    for measurement in measurements:
        f_pull, z_pull, f_release, z_release, title = ba.read_analyze(measurement, p)
        wlc, _ = func.WLC(f_pull, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])

        plt.ylabel('F (pN)')
        plt.xlabel('z (nm)')
        plt.tick_params(direction='in', top=True, right=True)
        plt.ylim(-1, 65)
        plt.xlim(500, 2300)

        plt.scatter(1000 * z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        plt.scatter(1000 * z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        plt.plot(wlc, f_pull, '--', color="black", label="WLC", zorder=10000)
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

    # TODO rotations
    # TODO calculate drift (for entire sample, use median)

    return


def main_fitfiles():

    fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\Fitfiles\\"

    fitfiles = []
    os.chdir(fitfile_path)
    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    for fitfile in fitfiles:

        f_pull, f_release, z_pull, z_release, z_fit_pull = ba.read_fitfiles(fitfile_path, fitfile, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)
        wlc, _ = func.WLC(f_wlc, L_bp=p['L_bp'], P_nm=p['P_nm'], S_pN=p['S_pN'])

        plt.ylabel('F (pN)')
        plt.xlabel('z (nm)')
        plt.tick_params(direction='in', top=True, right=True)

        plt.ylim(-1, 65)
        plt.xlim(0.5, 1.8)

        plt.scatter(z_pull, f_pull, color='darkgreen', label="Pull", s=10, zorder=25, facecolors='none')
        plt.scatter(z_release, f_release, color='lightgrey', s=10, zorder=15, label='Release', facecolors='none')
        plt.plot(wlc/1000, f_wlc, '--', color="black", label="WLC", zorder=100)
        plt.plot(z_fit_pull, f_pull, color='black', linewidth=3, label="Stat. Mech. Model fit", zorder=1000)

        plt.legend(loc=2, frameon=False)
        plt.title(fitfile[:-4]+" - fitfile")

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
    main_fitfiles()
    main_assemble_pars()
