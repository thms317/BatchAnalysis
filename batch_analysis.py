import numpy as np
from scipy.optimize import curve_fit
import functions as func
import matplotlib.pyplot as plt


def batch_analysis(beads, data_lines, headers, time, diff_magnet, force, title, file_name, analysis_path):
    # constants

    kBT = 4.114  # (pN nm) - Boltzmann factor
    Lc = 1500  # contour length (bp)
    p = 50  # persistence length (nm)
    S = 1000  # stretch modulus (pN)
    L = Lc * 0.34  # contour length (nm)
    x0 = 0  # offset (nm)

    for bead in range(0, beads):

        print("Processing bead " + str(bead) + " of " + str(beads))

        # load the data
        Z = []
        for x in data_lines:
            Z.append(float(x.split()[headers.index('Z' + str(bead) + ' (um)')]))

        # calculate drift for the individual bead
        slope = func.calc_drift_self(data_lines, headers, time, bead)

        # correcting drift
        Z_drift = []
        for n, t in enumerate(time):
            Z_drift.append(Z[n] - (slope / 1000) * t)
        Z = np.array(Z_drift)

        # split the data in pull/release-curve
        f_pull = []
        f_release = []
        z_pull = []
        z_release = []
        time_pull=[]
        time_release=[]

        trigger = []  # from what data point does the pulling trace start

        # if the differential of the magnet is positive -> pull, else -> release ('factor' since 0 does not work)
        for n, i in enumerate(diff_magnet):
            factor = max(diff_magnet / 1000)
            if i < -factor:
                trigger.append(n)
                f_pull.append(force[n])
                z_pull.append(Z[n])
                time_pull.append(time[n])
            if i > factor:
                f_release.append(force[n])
                z_release.append(Z[n])
                time_release.append(time[n])

        # wlc for reference
        wlc = []
        for f in f_pull:
            wlc.append(func.WLC(f, p, L, S, x0))

        # select data
        select_f = []
        select_z = []
        for n, f in enumerate(f_pull):
            if 20 < f < 30:
                select_f.append(f)
                select_z.append(Z[n + min(trigger)])

        # initial guesses
        x_init = 1

        # fit the WLC in fashion (x,y) - only fit offset, fix everything else
        popt, pcov = curve_fit(lambda f, x0: func.WLC(f, p, L, S, x0), select_f, select_z, p0=(x_init))
        std = np.sqrt(np.diag(pcov))  # returns the standard deviation

        x_fit = popt[0]

        z_pull -= x_fit  # subtract fitted offset from data
        z_release -= x_fit  # subtract fitted offset from data
        select_z -= x_fit

        a = np.percentile(z_pull, 1)
        dZ = "{0:.3f}".format(a - np.percentile(z_pull, 99))

        # plotting + saving

        # marker_size = 10
        #
        # plt.subplot(2, 1, 1)
        #
        # plt.title(str(title) + " / " + str(file_name) + " / bead " + str(bead) + " (dZ = " + str(dZ) + " nm)")
        # plt.scatter(time_pull, z_pull-a, facecolor='None', edgecolors="darkgreen", s=marker_size)
        # plt.scatter(time_release, z_release-a, facecolor='None', edgecolors="darkgrey", s=marker_size)
        # plt.ylim(0, 0.75)
        # plt.xlabel("Time (s)")
        # plt.ylabel("Extension ($\mu$m)")
        #
        # plt.subplot(2, 2, 3)
        #
        # plt.plot(wlc, f_pull, color='black', zorder=100)
        # plt.scatter(z_pull, f_pull, facecolor='None', edgecolors="darkgreen", s = marker_size)
        # plt.ylabel("Force (pN)")
        # plt.xlabel("Extension ($\mu$m)")
        # plt.xlim(0, 0.75)
        #
        # plt.subplot(2, 2, 4)
        #
        # # plt.plot(wlc, f_pull, color='black', zorder=100)
        # plt.scatter(z_release, f_release, facecolor='None', edgecolors="darkgrey", s = marker_size)
        # plt.ylabel("Force (pN)")
        # plt.xlabel("Extension ($\mu$m)")
        # plt.xlim(0, 0.75)
        #
        # plt.savefig(analysis_path + "dZ_" + str(dZ) + "_" +  str(title) + "_" + str(file_name) + "_bead" + str(bead) + '_subplot.png', dpi=300, bbox_inches='tight')
        #
        # plt.close()

    return
