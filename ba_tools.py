import functions as func
import numpy as np
from scipy.optimize import curve_fit
import os
import pandas as pd


def read_analyze(measurement, pars):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    global_drift = p['global_drift']

    # constants from parameters
    L_bp = p['L_bp']  # contour length (bp)
    S_pN = p['S_pN']  # stretch modulus (pN)
    P_nm = p['P_nm']  # persistence length (nm)

    #  laptop or work PC
    folder = "C:\\Users\\tbrouwer\\Desktop\\Data\\"
    bool = os.path.isdir(folder)
    if bool == False:
        folder = "C:\\Users\\brouw\\Desktop\\Data\\"

    # open data
    file_location = folder + str(measurement[0]) + "\\"
    file_name = "data_" + str(measurement[1])
    file_extension = ".dat"
    file_all = file_location + file_name + file_extension
    bead = int(measurement[2])

    # title
    title = int(func.get_num(file_location))

    # read DataFrame
    df = pd.read_csv(file_all, sep="\t")

    magnet = np.array(df['Stepper shift (mm)'])
    time = np.array(df['Time (s)'])
    force = func.calc_force(magnet)
    Z = np.array(df['Z' + str(bead) + ' (um)'])

    if global_drift:

        # number of beads
        file_log = file_location + file_name + ".log"
        f = open(file_log, 'r')
        beads = f.readlines()[9]
        f.close()

        beads = int(func.get_num(beads))

        drift = []
        for i in range(beads):
            z_drift = np.array(df['Z' + str(i) + ' (um)'])
            amplitude_drift = np.array(df['Amp' + str(i) + ' (a.u.)'])

            # does the tether rupture?
            rupt = func.rupture(time, amplitude_drift)

            if rupt == False:
                drift.append(func.drift_self(z_drift, time))

        drift_med = np.median(drift)
        Z = Z - (drift_med / 1000) * time
        p['drift'] = str(round(drift_med, 3)) + " nm/s (global drift)"

    else:

        amplitude = np.array(df['Amp' + str(bead) + ' (a.u.)'])

        # does the tether rupture?
        rupt = func.rupture(time, amplitude)

        if rupt == False:
            # correcting the drift using self
            drift = func.drift_self(Z, time)
            Z = Z - (drift / 1000) * time
            p['drift'] = str(round(drift, 3)) + " nm/s"
        else:
            p['drift'] = 'uncorrected'

    # calculating the first derivative of magnet
    dx = np.diff(time)
    dy = np.diff(magnet)
    diff_magnet = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # split in pull & release
    factor = max(diff_magnet / 1000)
    time_pull = time[diff_magnet < factor]
    f_pull = force[diff_magnet < factor]
    f_release = force[diff_magnet > factor]
    z_pull = Z[diff_magnet < factor]
    z_release = Z[diff_magnet > factor]

    # select high-force data
    select_f = f_pull[np.where((f_pull > 30) & (f_pull < 50) & (time_pull > 60) & (time_pull < 140))]
    select_z = z_pull[np.where((f_pull > 30) & (f_pull < 50) & (time_pull > 60) & (time_pull < 140))]

    # fit the WLC in fashion (x,y) - only fit offset, fix everything else
    popt, pcov = curve_fit(lambda f, z0: func.WLC_fit(f, P_nm, L_bp * 0.34, S_pN, z0), select_f, select_z, p0=1)

    z_fit = popt[0]

    # subtract fitted offset from data
    z_pull -= z_fit
    z_release -= z_fit
    select_z -= z_fit

    title = str(title) + '_' + str(measurement[1]) + '_' + str(measurement[2]) + '_' + str(measurement[3])

    return f_pull, z_pull, f_release, z_release, title


def read_analyze_rot(measurement, pars):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    global_drift = p['global_drift']

    # constants from parameters
    L_bp = p['L_bp']  # contour length (bp)
    S_pN = p['S_pN']  # stretch modulus (pN)
    P_nm = p['P_nm']  # persistence length (nm)

    #  laptop or work PC
    folder = "C:\\Users\\tbrouwer\\Desktop\\Data\\"
    bool = os.path.isdir(folder)
    if bool == False:
        folder = "C:\\Users\\brouw\\Desktop\\Data\\"

    # open data
    file_location = folder + str(measurement[0]) + "\\"
    file_name = "data_" + str("{0:0=3d}".format(int(measurement[1]) - 1))
    file_extension = ".dat"
    file_all = file_location + file_name + file_extension
    bead = int(measurement[2])

    # title
    title = int(func.get_num(file_location))

    # read DataFrame
    df = pd.read_csv(file_all, sep="\t")

    try:
        rotation = np.array(df['Stepper rot (turns)'])
    except:
        print("No rotation measurement")
        return [], [], [], [], []

    time = np.array(df['Time (s)'])
    Z = np.array(df['Z' + str(bead) + ' (um)'])

    if global_drift:

        # number of beads
        file_log = file_location + file_name + ".log"
        f = open(file_log, 'r')
        beads = f.readlines()[9]
        f.close()

        beads = int(func.get_num(beads))

        drift = []
        for i in range(beads):
            z_drift = np.array(df['Z' + str(i) + ' (um)'])
            amplitude_drift = np.array(df['Amp' + str(i) + ' (a.u.)'])

            # does the tether rupture?
            rupt = func.rupture(time, amplitude_drift)

            if rupt == False:
                drift.append(func.drift_self(z_drift, time))

        drift_med = np.median(drift)
        Z = Z - (drift_med / 1000) * time
        p['drift'] = str(round(drift_med, 3)) + " nm/s (global drift)"

    else:

        amplitude = np.array(df['Amp' + str(bead) + ' (a.u.)'])

        # does the tether rupture?
        rupt = func.rupture(time, amplitude)

        if rupt == False:
            # correcting the drift using self
            drift = func.drift_self(Z, time)
            Z = Z - (drift / 1000) * time
            p['drift'] = str(round(drift, 3)) + " nm/s"
        else:
            p['drift'] = 'uncorrected'

    # calculating the first derivative of magnet
    dx = np.diff(time)
    dy = np.diff(rotation)
    diff_rot = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # offset the Z
    Z -= np.mean(Z)

    # calculate LND
    LND = rotation / (p['L_bp'] / 10.4)

    # split in pull & release
    factor = max(diff_rot / 1000)
    twist_pos = rotation[diff_rot < factor]
    twist_neg = rotation[diff_rot > factor]
    z_pos = Z[diff_rot < factor]
    z_neg = Z[diff_rot > factor]
    lnd_pos = LND[diff_rot < factor]
    lnd_neg = LND[diff_rot > factor]

    return twist_pos, twist_neg, z_pos, z_neg, lnd_pos, lnd_neg

def read_fitfiles(fitfile_path, fitfile, pars):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    mask = False
    standard_trajectory = False

    #  laptop or work PC
    bool = os.path.isdir(fitfile_path)
    if bool == False:
        fitfile_path = fitfile_path.replace("tbrouwer", "brouw")

    # open data
    file_all = fitfile_path + fitfile

    # read DataFrame
    df = pd.read_csv(file_all, sep="\t")

    time = np.array(df['t (s)'])
    force = np.array(df['F (pN)'])
    z = np.array(df['z (um)'])
    z_fit = np.array(df['z fit (um)'])

    # transitions
    trans_number = int(df.columns[-1][3:])  # number of transitions
    t1 = df.columns.get_loc("T1_0")  # locations
    t2 = df.columns.get_loc("T2_0")
    t3 = df.columns.get_loc("T3_0")
    T1 = np.array(df.iloc[:, t1:1 + trans_number + t1])  # transitions
    T2 = np.array(df.iloc[:, t2:1 + trans_number + t2])
    T3 = np.array(df.iloc[:, t3:1 + trans_number + t3])

    if mask:
        # calculating the mask
        # rupt, mask, test = filter_rupture(fitfile, test=True)  # use to check if results do not make sense
        rupt, mask = filter_rupture(fitfile)

        # applying the mask
        time = time[mask == 1]
        force = force[mask == 1]
        z = z[mask == 1]
        z_fit = z_fit[mask == 1]

        T1 = T1[mask == 1]  # filter the transitions
        T2 = T2[mask == 1]
        T3 = T3[mask == 1]

    transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

    # calculating the first derivative of force
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # split in pull & release
    factor = max(diff_force / 1000)

    if standard_trajectory:
        f_pull = force[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        f_release = force[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        z_pull = z[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        z_release = z[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 50) & (time < 80))]

    else:
        f_pull = force[np.where(diff_force > factor)]
        f_release = force[np.where(diff_force < factor)]
        z_pull = z[np.where(diff_force > factor)]
        z_release = z[np.where(diff_force < factor)]
        z_fit_pull = z_fit[np.where(diff_force > factor)]

    return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions, force


def filter_rupture(fitfile, test=False):
    #  laptop or work PC
    folder = "C:\\Users\\tbrouwer\\Desktop\\Data\\"
    bool = os.path.isdir(folder)
    if bool == False:
        folder = "C:\\Users\\brouw\\Desktop\\Data\\"

    # open data
    file_location = folder + str(fitfile[:6]) + "\\"
    file_name = str(fitfile[7:15]) + ".dat"
    bead = fitfile[16:-4]

    # because of this annoying bug in the LabVIEW program, offset bead no. by 2
    bead = int(bead)
    bead = str(bead + 2)

    # title
    title = fitfile[:-4]

    # read DataFrame
    df = pd.read_csv(file_location + file_name, sep="\t")

    time = np.array(df['Time (s)'])
    amplitude = np.array(df['Amp' + str(bead) + ' (a.u.)'])
    Z = np.array(df['Z' + str(bead) + ' (um)'])

    # does the tether rupture?
    rupt, mask = func.rupture(time, amplitude, mask=True)

    if test:
        return rupt, mask, Z
    else:
        return rupt, mask


def build_measurements(file_location, p):
    xl = pd.ExcelFile(file_location)
    df = xl.parse("Sheet1")
    selected = np.array(df['Selected'])

    data = np.array(df['File'])[np.where(selected == 1)]
    bead = np.array(df['Trace'])[np.where(selected == 1)]
    date = np.full((len(bead)), file_location[-11:-5])
    type = np.full((len(bead)), p['NRL'])

    measurements = []
    for n in range(len(bead)):
        measurements.append([date[n], data[n][5:8], bead[n], type[n]])

    return measurements
