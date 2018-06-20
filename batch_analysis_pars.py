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
    pars['L_bp'] = 4753  # number of base pairs
    pars['P_nm'] = 50  # persistence length
    pars['S_pN'] = 1000  # stretch modulus
    pars['z0_nm'] = 0  # offset in nm / subunit
    pars['NRL'] = 168  # nucleosome repeat length
    pars['repeats'] = 16  # number of repeats
    pars['type'] = "Human"  # type of histone
    pars['NRL_str'] = str(pars['NRL'])+'x'+str(pars['repeats'])+'_'+pars['type']  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    return pars

p = default_pars()

def main_pars():
    fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\"

    fitfiles = []
    os.chdir(fitfile_path)
    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    ass_fit_pars = []
    ass_fit_errors = []

    for fitfile in fitfiles:
        print("Processing fitfile... " + str(fitfile))

        # read pars from logfile
        logfile = fitfile[:-3] + "log"
        fit_pars, fit_errors, table = ba.read_logfile(fitfile_path, logfile)
        ass_fit_pars.append(fit_pars)
        ass_fit_errors.append(fit_errors)

    ass_fit_pars=np.transpose(ass_fit_pars)

    plt.scatter(ass_fit_pars[0],ass_fit_pars[3])
    plt.xlabel("Number of Nucleosomes")
    plt.ylabel("G1 (kT)")
    plt.show()

    return

if __name__ == "__main__":
    main_pars()
