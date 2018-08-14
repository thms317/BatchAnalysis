import os
import matplotlib.pyplot as plt
import functions as func
import ba_tools as ba
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
    pars['NRL'] = 170  # nucleosome repeat length
    pars['repeats'] = 16  # number of repeats
    pars['type'] = "human"  # type of histone
    pars['NRL_str'] = str(pars['NRL'])+'x'+str(pars['repeats'])+'_'+pars['type']  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    pars['standard'] = True
    pars['radius_um'] = []  # radius of circle (um)
    return pars

p = default_pars()


plt.close("all")

# table_path = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Selection Tables\\"
table_path = "C:\\Users\\tbrouwer\\Desktop\\Data\\"
table_file = "180814_170"

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