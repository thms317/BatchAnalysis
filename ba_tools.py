import functions as func
import numpy as np
from scipy.optimize import curve_fit
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def read_analyze(measurement, pars, data_path):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    global_drift = True

    # constants from parameters
    L_bp = p['L_bp']  # contour length (bp)
    S_pN = p['S_pN']  # stretch modulus (pN)
    P_nm = p['P_nm']  # persistence length (nm)

    # open data
    file_location = data_path
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

    # number of beads
    file_log = file_location + file_name + ".log"
    f = open(file_log, 'r')
    try:
        beads = f.readlines()[9]
    except:
        headers = list(df)
        # get number of beads
        beads = headers[len(headers) - 1]
    f.close()

    try:
        beads = int(func.get_num(beads))
    except:
        beads = func.get_int(beads)

    drift = []
    for i in range(beads):
        z_drift = np.array(df['Z' + str(i) + ' (um)'])
        try:
            amplitude_drift = np.array(df['Amp' + str(i) + ' (a.u.)'])
        except:
            amplitude_drift = np.array(df['Amp a.u.'])

        # does the tether rupture?
        rupt, _ = func.rupture(time, amplitude_drift)

        if rupt == False:
            drift.append(func.drift_self(z_drift, time))

    if not drift:
        p['drift'] = 'uncorrected'
        drift_med = 0
    else:
        drift_med = float(np.median(drift))
        Z = Z - (drift_med / 1000) * time
        p['drift'] = str(round(drift_med, 3)) + " nm/s (global drift)"

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
    select_f = f_pull[np.where((f_pull > 50) & (f_pull < 55) & (time_pull > 60) & (time_pull < 140))]
    select_z = z_pull[np.where((f_pull > 50) & (f_pull < 55) & (time_pull > 60) & (time_pull < 140))]

    if select_f.size != 0:

        # fit the WLC in fashion (x,y) - only fit offset, fix everything else
        popt, pcov = curve_fit(lambda f, z0: func.WLC_fit(f, P_nm, L_bp * 0.34, S_pN, z0), select_f, select_z, p0=1)

        z_fit = popt[0]

        # subtract fitted offset from data
        z_pull -= z_fit
        z_release -= z_fit
        select_z -= z_fit

    else:
        z_fit = 0

    title = str(title) + '_' + str(measurement[1]) + '_' + str(measurement[2]) + '_' + str(measurement[3])

    return f_pull, z_pull, f_release, z_release, title, drift_med, z_fit


def read_analyze_rot(measurement, pars, data_path):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    # constants from parameters
    L_bp = p['L_bp']  # contour length (bp)
    S_pN = p['S_pN']  # stretch modulus (pN)
    P_nm = p['P_nm']  # persistence length (nm)

        # open data
    file_location = data_path
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
        p['drift'] = 'uncorrected'
        return [], [], [], [], [], []

    time = np.array(df['Time (s)'])
    X = np.array(df['X' + str(bead) + ' (um)'])
    Y = np.array(df['Y' + str(bead) + ' (um)'])
    Z = np.array(df['Z' + str(bead) + ' (um)'])

    # number of beads
    file_log = file_location + file_name + ".log"
    f = open(file_log, 'r')
    beads = f.readlines()[9]
    f.close()

    beads = int(func.get_num(beads))

    drift = []
    for i in range(beads):
        z_drift = np.array(df['Z' + str(i) + ' (um)'])
        try:
            amplitude_drift = np.array(df['Amp' + str(i) + ' (a.u.)'])
        except:
            amplitude_drift = np.array(df['Amp a.u.'])

        # does the tether rupture?
        rupt, _ = func.rupture(time, amplitude_drift)

        if rupt == False:
            drift.append(func.drift_self(z_drift, time))

    if not drift:
        p['drift'] = 'uncorrected'

    else:
        drift_med = float(np.median(drift))
        Z = Z - (drift_med / 1000) * time
        p['drift'] = str(round(drift_med, 3)) + " nm/s (global drift)"

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
    x_pos = X[diff_rot < factor]
    x_neg = X[diff_rot > factor]
    y_pos = Y[diff_rot < factor]
    y_neg = Y[diff_rot > factor]
    z_pos = Z[diff_rot < factor]
    z_neg = Z[diff_rot > factor]
    lnd_pos = LND[diff_rot < factor]
    lnd_neg = LND[diff_rot > factor]

    x_sel = np.concatenate((x_pos,x_neg))
    y_sel = np.concatenate((y_pos,y_neg))

    # fit circle to extract radius
    x_m, y_m, R = fit_circle(x_sel,y_sel)
    # alternative method (currently works better)
    R2 = (abs(np.percentile(x_sel,1)-np.percentile(x_sel,99))/2+abs(np.percentile(y_sel,1)-np.percentile(y_sel,99))/2)/2
    x_m_2 = np.percentile(x_sel,1)+R2
    y_m_2 = np.percentile(y_sel,1)+R2

    p['radius_um'] = R2

    # plt.close()
    # plt.scatter(x_sel,y_sel)
    # plt.scatter(x_m, y_m)
    # plt.scatter(x_m_2, y_m_2, color='r')
    # plt.plot(x_m+R,y_m)
    # plt.show()

    return twist_pos, twist_neg, z_pos, z_neg, lnd_pos, lnd_neg


def read_fitfiles(fitfile_path, fitfile, pars, meas_pars):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    standard_trajectory = p['standard']

    m = meas_pars

    m['date'] = fitfile[:6]
    m['data'] = fitfile[12:15]
    m['bead'] = fitfile[:-4][16:]
    m['repeats'] = p['repeats']

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

    # calculating the mask
    rupt, mask = filter_rupture(fitfile, m)
    m['points'] = len(time)
    m['points_frac'] = m['points'] / m['points_exp']

    # applying the mask
    time = time[mask == 1]
    force = force[mask == 1]
    z = z[mask == 1]
    z_fit = z_fit[mask == 1]

    T1 = T1[mask == 1]  # filter the transitions
    T2 = T2[mask == 1]
    T3 = T3[mask == 1]

    m['dZ_um'] = abs(np.percentile(z, 1)-np.percentile(z, 99))

    # calculating the first derivative of force
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # split in pull & release
    factor = max(diff_force / 1000)

    if standard_trajectory:

        f_pull = force[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        f_release = force[np.where((diff_force < factor) & (time > 120))]
        z_pull = z[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        z_release = z[np.where((diff_force < factor) & (time > 120))]
        time_pull = time[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        time_release = time[np.where((diff_force < factor) & (time > 120))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T1 = T1[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T2 = T2[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T3 = T3[np.where((diff_force > factor) & (time > 90) & (time < 130))]

        # f_pull = force[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # f_release = force[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        # z_pull = z[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # z_release = z[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        # time_pull = time[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # time_release = time[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        # z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # T1 = T1[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # T2 = T2[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        # T3 = T3[np.where((diff_force > factor) & (time > 50) & (time < 80))]

    else:
        f_pull = force[np.where(diff_force > factor)]
        f_release = force[np.where(diff_force < factor)]
        z_pull = z[np.where(diff_force > factor)]
        z_release = z[np.where(diff_force < factor)]
        z_fit_pull = z_fit[np.where(diff_force > factor)]
        T1 = T1[np.where(diff_force > factor)]
        T2 = T2[np.where(diff_force > factor)]
        T3 = T3[np.where(diff_force > factor)]

    transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

    return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions
    # return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions, time_pull, time_release, force, z, time

def read_fitfiles_plain(fitfile_path, fitfile, pars, standard_trajectory=False, evaluate_ruptures=False):

    p = pars

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
    # z = np.array(df['selected z (um)'])
    z_fit = np.array(df['z fit (um)'])

    # transitions
    trans_number = int(df.columns[-1][3:])  # number of transitions
    t1 = df.columns.get_loc("T1_0")  # locations
    t2 = df.columns.get_loc("T2_0")
    t3 = df.columns.get_loc("T3_0")
    T1 = np.array(df.iloc[:, t1:1 + trans_number + t1])  # transitions
    T2 = np.array(df.iloc[:, t2:1 + trans_number + t2])
    T3 = np.array(df.iloc[:, t3:1 + trans_number + t3])

    if evaluate_ruptures:

        # calculating the mask
        rupture, rupt_index, mask = filter_rupture_Z(time, force, z)

        # applying the mask
        time = time[mask == 1]
        force = force[mask == 1]
        z = z[mask == 1]
        z_fit = z_fit[mask == 1]

        T1 = T1[mask == 1]  # filter the transitions
        T2 = T2[mask == 1]
        T3 = T3[mask == 1]

    # calculating the first derivative of force
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # split in pull & release
    factor = max(diff_force / 1000)

    if p['NRL'] == 168 or p['NRL'] == 169:
    # if p['NRL'] == 169:
        f_pull = force[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        f_release = force[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        z_pull = z[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        z_release = z[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        time_pull = time[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        time_release = time[np.where((diff_force < factor) & (time > 75) & (time < 125))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        T1 = T1[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        T2 = T2[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        T3 = T3[np.where((diff_force > factor) & (time > 50) & (time < 80))]

        transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

        return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions

    if standard_trajectory:

        f_pull = force[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        f_release = force[np.where((diff_force < factor) & (time > 120))]
        z_pull = z[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        z_release = z[np.where((diff_force < factor) & (time > 120))]
        time_pull = time[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        time_release = time[np.where((diff_force < factor) & (time > 120))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T1 = T1[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T2 = T2[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T3 = T3[np.where((diff_force > factor) & (time > 90) & (time < 130))]

    else:
        f_pull = force[np.where(diff_force > factor)]
        f_release = force[np.where(diff_force < factor)]
        z_pull = z[np.where(diff_force > factor)]
        z_release = z[np.where(diff_force < factor)]
        z_fit_pull = z_fit[np.where(diff_force > factor)]
        T1 = T1[np.where(diff_force > factor)]
        T2 = T2[np.where(diff_force > factor)]
        T3 = T3[np.where(diff_force > factor)]

    transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

    return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions

def read_fitfiles_assembly(fitfile_path, fitfile, pars, standard_trajectory=False, evaluate_ruptures=False):

    p = pars

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
    # z = np.array(df['z (um)'])
    z = np.array(df['selected z (um)'])
    z_fit = np.array(df['z fit (um)'])

    # transitions
    trans_number = int(df.columns[-1][3:])  # number of transitions
    t1 = df.columns.get_loc("T1_0")  # locations
    t2 = df.columns.get_loc("T2_0")
    t3 = df.columns.get_loc("T3_0")
    T1 = np.array(df.iloc[:, t1:1 + trans_number + t1])  # transitions
    T2 = np.array(df.iloc[:, t2:1 + trans_number + t2])
    T3 = np.array(df.iloc[:, t3:1 + trans_number + t3])

    if evaluate_ruptures:

        # calculating the mask
        rupture, rupt_index, mask = filter_rupture_Z(time, force, z)

        # applying the mask
        time = time[mask == 1]
        force = force[mask == 1]
        z = z[mask == 1]
        z_fit = z_fit[mask == 1]

        T1 = T1[mask == 1]  # filter the transitions
        T2 = T2[mask == 1]
        T3 = T3[mask == 1]

    # calculating the first derivative of force
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element

    # split in pull & release
    factor = max(diff_force / 1000)

    if standard_trajectory:

        f_pull = force[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        f_release = force[np.where((diff_force < factor) & (time > 120))]
        z_pull = z[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        z_release = z[np.where((diff_force < factor) & (time > 120))]
        time_pull = time[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        time_release = time[np.where((diff_force < factor) & (time > 120))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T1 = T1[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T2 = T2[np.where((diff_force > factor) & (time > 90) & (time < 130))]
        T3 = T3[np.where((diff_force > factor) & (time > 90) & (time < 130))]

    else:
        f_pull = force[np.where(z != np.nan)]
        f_release = force[np.where(z != np.nan)]
        z_pull = z[np.where(z != np.nan)]
        z_release = z[np.where(z != np.nan)]
        z_fit_pull = z_fit[np.where(z != np.nan)]
        T1 = T1[np.where(z != np.nan)]
        T2 = T2[np.where(z != np.nan)]
        T3 = T3[np.where(z != np.nan)]

    transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

    return f_pull, f_release, z_pull, z_release, z_fit_pull, transitions


def filter_rupture(fitfile, meas_pars, test=False):
    m = meas_pars

    #  laptop or work PC
    folder = "C:\\Users\\tbrouwer\\Desktop\\Data\\"
    bool = os.path.isdir(folder)
    if bool == False:
        folder = "C:\\Users\\brouw\\Desktop\\Data\\"

    # open data
    file_location = folder + str(fitfile[:6]) + "\\"
    file_name = str(fitfile[7:15]) + ".dat"
    bead = fitfile[16:-4]

    # read logfile
    log_name = str(fitfile[7:15]) + ".log"
    f = open(file_location + log_name, 'r')
    log = f.readlines()[:]
    f.close()
    m['beads'] = func.get_num(log[9])
    framerate = func.get_num(log[5])

    # because of this annoying bug in the LabVIEW program, offset bead no. by 2
    bead = int(bead)
    # bead = str(bead + 2)
    bead = str(bead)

    # title
    title = fitfile[:-4]

    # read DataFrame
    df = pd.read_csv(file_location + file_name, sep="\t")

    time = np.array(df['Time (s)'])
    amplitude = np.array(df['Amp' + str(bead) + ' (a.u.)'])
    X = np.array(df['X' + str(bead) + ' (um)'])
    Y = np.array(df['Y' + str(bead) + ' (um)'])
    Z = np.array(df['Z' + str(bead) + ' (um)'])

    m['X0_um'] = np.mean(X[np.logical_not(np.isnan(X))])
    m['Y0_um'] = np.mean(Y[np.logical_not(np.isnan(Y))])
    m['Z0_um'] = np.mean(Z[np.logical_not(np.isnan(Z))])
    m['amp'] = np.mean(amplitude[np.logical_not(np.isnan(amplitude))])
    m['points_exp'] = int(time[-1]*framerate)

    # does the tether rupture?
    rupt, peak_index, mask = func.rupture(time, amplitude, mask=True)

    if rupt:
        magnet = np.array(df['Stepper shift (mm)'])
        force = func.calc_force(magnet)
        f_rupt = force[peak_index]
        m['f_rupt_pN'] = f_rupt

    if test:
        return rupt, mask, Z
    else:
        return rupt, mask

def filter_rupture_Z(time, force, Z):

    # calculating the first derivative of amplitude
    dx = np.diff(time)
    dy = np.diff(Z)
    diff_z = np.append([0], np.divide(dy, dx))
    diff_z -= np.mean(diff_z)

    dx = np.diff(time)
    dy = np.diff(diff_z)
    diff2_z = abs(np.append([0], np.divide(dy, dx)))
    diff2_z -= np.mean(diff_z)

    rupture = False

    for n, ddz in enumerate(diff2_z):
        rupt_index = n
        if ddz > 10 * np.median(diff2_z):
            rupture = True
            break

    mask = np.ones(len(Z))

    if rupture:
        mask_on = np.ones(rupt_index)
        fuck_it_mask_off = np.zeros(len(Z) - rupt_index)
        mask = np.concatenate((mask_on, fuck_it_mask_off))

    return rupture, rupt_index, mask


def build_measurements(folder, datx, pars):
    try:
        p = pars
    except:
        print('Error: no parameters')
        return

    # read DataFrame
    df = pd.read_csv(folder + datx, sep="\t")

    selected = np.array(df['Selected'])

    data = np.array(df['File'])[np.where(selected == 1)]
    bead = np.array(df['Trace'])[np.where(selected == 1)]
    date, type = [], []
    for n in range(len(bead)):
        date.append(datx[:6])
        type.append(p['NRL_str'])

    measurements = []
    for n in range(len(bead)):
        measurements.append([date[n], data[n][5:8], bead[n], type[n]])

    return measurements


def read_logfile(logfile_path, logfile, meas_pars):

    m = meas_pars

    f = open(logfile_path + logfile, 'r')
    log = f.readlines()[:]
    f.close()

    for n, line in enumerate(log):
        log[n] = line.rstrip()

    i = log.index("[Fit parameters]")

    fit_pars = []
    fit_pars.append(func.get_num(log[i + 7]))  # N_nuc
    fit_pars.append(func.get_num(log[i + 10][14:]))  # N_unfolded
    fit_pars.append(round(func.get_num(log[i + 9]),3))  # k
    fit_pars.append(round(func.get_num(log[i + 13][4:]),3))  # G1
    fit_pars.append(round(func.get_num(log[i + 14][4:]),3))  # G2
    fit_pars.append(func.get_num(log[i + 16][16:]))  # degeneracy
    fit_pars.append(func.get_num(log[i + 7]) - func.get_num(log[i + 10][14:])) # stacked nucleosomes
    fit_pars.append(func.get_num(log[i + 5])) # Stretch Modulus

    m['drift_nm_s'] = func.get_num(log[i + 2])
    m['L_bp'] = func.get_num(log[i + 3])
    m['P_nm'] = func.get_num(log[i + 4])
    m['S_pN'] = func.get_num(log[i + 5])
    m['NRL'] = func.get_num(log[i + 6])
    m['N_nuc'] = func.get_num(log[i + 7])
    m['N_tet'] = func.get_num(log[i + 10][14:])
    m['N_stack'] = func.get_num(log[i + 7]) - func.get_num(log[i + 10][14:])
    m['k_pN_nm'] = func.get_num(log[i + 9])
    m['G1_kT'] = func.get_num(log[i + 13][4:])
    m['G2_kT'] = func.get_num(log[i + 14][4:])
    m['G3_kT'] = func.get_num(log[i + 15][4:])
    m['L_fold'] = func.get_num(log[i + 8])
    m['L_unwrap'] = func.get_num(log[i + 11])
    m['L_ext'] = func.get_num(log[i + 12])
    m['degeneracy'] = func.get_num(log[i + 16][16:])

    errors = [[], [], []]
    try:
        j = log.index("[Fit local error]")

        for error in range(j, len(log)):
            if "k folded (pN/nm)" in log[error]:
                errors[0] = round(func.get_num(log[error]),3)  # place k-error in errors[0]
                m['k_pN_nm_se'] = func.get_num(log[error])
            if "G1 (kT)" in log[error]:
                strip_log = log[error][2:]
                errors[1] = round(func.get_num(strip_log),3)  # place G1-error in errors[1]
                m['G1_kT_se'] = func.get_num(strip_log)
            if "G2 (kT)" in log[error]:
                strip_log = log[error][2:]
                errors[2] = round(func.get_num(strip_log),3)  # place G2-error in errors[2]
                m['G2_kT_se'] = func.get_num(strip_log)

        p0 = "N_nuc = " + str(fit_pars[0])
        p1 = "N_unfolded = " + str(fit_pars[1])
        p2 = "k (pN/nm) = " + str(fit_pars[2]) + " +/- " + str(errors[0])
        p3 = "G1 (kT) = " + str(fit_pars[3]) + " +/- " + str(errors[1])
        p4 = "G2 (kT) = " + str(fit_pars[4]) + " +/- " + str(errors[2])
        p5 = "Degeneracy = " + str(fit_pars[5])
        p6 = "N_stacked = " + str(fit_pars[6])
        if int(fit_pars[7]) != 1000:
            p7 = "S (pN) = " + str(fit_pars[7])
        else:
            p7 = ''

    except:

        p0 = "N_nuc = " + str(fit_pars[0])
        p1 = "N_unfolded = " + str(fit_pars[1])
        p2 = "k (pN/nm) = " + str(fit_pars[2])
        p3 = "G1 (kT) = " + str(fit_pars[3])
        p4 = "G2 (kT) = " + str(fit_pars[4])
        p5 = "Degeneracy = " + str(fit_pars[5])
        p6 = "Stacked Nucleosomes = " + str(fit_pars[6])
        if int(fit_pars[7]) != 1000:
            p7 = "S (pN) = " + str(fit_pars[7])
        else:
            p7 = ''

    table = [p0, p1, p2, p3, p4, p5, p6, p7]

    return fit_pars, errors, table

def read_logfile_clean(logfile_path, logfile, pars):

    p = pars

    p['date'] = logfile[:6]
    p['data'] = logfile[12:15]
    p['bead'] = logfile[:-4][16:]

    f = open(logfile_path + logfile, 'r')
    log = f.readlines()[:]
    f.close()

    for n, line in enumerate(log):
        log[n] = line.rstrip()

    i = log.index("[Fit parameters]")

    fit_pars = []
    fit_pars.append(func.get_num(log[i + 7]))  # N_nuc
    fit_pars.append(func.get_num(log[i + 10][14:]))  # N_unfolded
    fit_pars.append(round(func.get_num(log[i + 9]),3))  # k
    fit_pars.append(round(func.get_num(log[i + 13][4:]),3))  # G1
    fit_pars.append(round(func.get_num(log[i + 14][4:]),3))  # G2
    fit_pars.append(func.get_num(log[i + 16][16:]))  # degeneracy
    fit_pars.append(func.get_num(log[i + 7]) - func.get_num(log[i + 10][14:])) # stacked nucleosomes
    fit_pars.append(func.get_num(log[i + 5])) # Stretch Modulus

    p['drift_nm_s'] = func.get_num(log[i + 2])
    p['L_bp'] = func.get_num(log[i + 3])
    p['P_nm'] = func.get_num(log[i + 4])
    p['S_pN'] = func.get_num(log[i + 5])
    p['NRL'] = func.get_num(log[i + 6])
    p['N_nuc'] = func.get_num(log[i + 7])
    p['N_tet'] = func.get_num(log[i + 10][14:])
    p['N_stack'] = func.get_num(log[i + 7]) - func.get_num(log[i + 10][14:])
    p['k_pN_nm'] = func.get_num(log[i + 9])
    p['G1_kT'] = func.get_num(log[i + 13][4:])
    p['G2_kT'] = func.get_num(log[i + 14][4:])
    p['G3_kT'] = func.get_num(log[i + 15][4:])
    p['L_fold'] = func.get_num(log[i + 8])
    p['L_unwrap'] = func.get_num(log[i + 11])
    p['L_ext'] = func.get_num(log[i + 12])
    p['degeneracy'] = func.get_num(log[i + 16][16:])

    errors = [[], [], []]
    try:
        j = log.index("[Fit local error]")

        for error in range(j, len(log)):
            if "k folded (pN/nm)" in log[error]:
                errors[0] = round(func.get_num(log[error]),3)  # place k-error in errors[0]
                p['k_pN_nm_se'] = func.get_num(log[error])
            if "G1 (kT)" in log[error]:
                strip_log = log[error][2:]
                errors[1] = round(func.get_num(strip_log),3)  # place G1-error in errors[1]
                p['G1_kT_se'] = func.get_num(strip_log)
            if "G2 (kT)" in log[error]:
                strip_log = log[error][2:]
                errors[2] = round(func.get_num(strip_log),3)  # place G2-error in errors[2]
                p['G2_kT_se'] = func.get_num(strip_log)

    except:
        pass

    return fit_pars, errors


def plot_hist(ass_fit_pars, ass_fit_errors, new_path, p, show_plot = True):

    ass_fit_pars = np.transpose(np.array(ass_fit_pars))
    # ass_fit_errors = np.transpose(np.array(ass_fit_errors))

    fig = plt.figure(figsize=(30, 18))
    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', linewidth=3)

    # # number of nucleosomes
    # ax0 = fig.add_subplot(2, 2, 1)
    #
    # ax0.set_ylabel('count')
    # ax0.set_xlabel('N_nuc')
    # ax0.tick_params(direction='in', top=True, right=True)
    # ax0.set_title("Number of Nucleosomes")
    # ax0.hist(ass_fit_pars[0])
    #
    # # number of tetrasomes
    # ax1 = fig.add_subplot(2, 2, 2)
    #
    # ax1.set_ylabel('count')
    # ax1.set_xlabel('N_unfolded')
    # ax1.tick_params(direction='in', top=True, right=True)
    # ax1.set_title("Number of Tetrasomes")
    # ax1.hist(ass_fit_pars[1])

    # # number of nucleosomes
    # ax0 = fig.add_subplot(2, 2, 1)
    #
    # ax0.set_ylabel('count')
    # ax0.set_xlabel('N_nuc')
    # ax0.tick_params(direction='in', top=True, right=True)
    # ax0.set_title("Number of Nucleosomes")
    # ax0.hist(ass_fit_pars[0])

    # number of stacked nucleosomes in fiber
    ax1 = fig.add_subplot(2, 2, 1)

    ax1.set_ylabel('count')
    ax1.set_xlabel('Stacked Nucleosomes in Fiber')
    ax1.tick_params(direction='in', top=True, right=True)
    ax1.set_title("Stacked Nucleosomes")
    stacked=[]
    for n,x in enumerate(ass_fit_pars[0]):
        stacked.append(ass_fit_pars[0][n]-ass_fit_pars[1][n])
    binwidth = 1 # nucleosomes
    ax1.hist(stacked, bins=np.arange(min(stacked), max(stacked) + binwidth, binwidth), edgecolor='black', linewidth=1.2)
    ax1.set_xlim(0,30)
    ax1.set_ylim(0,25)

    # stiffness

    ax2 = fig.add_subplot(2, 2, 3)

    ax2.set_ylabel('count')
    ax2.set_xlabel('k (pN/nm)')
    ax2.tick_params(direction='in', top=True, right=True)
    ax2.set_title("Fiber Stiffness")
    # binwidth = 0.05  # pN/nm - one start
    binwidth = 0.1  # pN/nm - two start
    ax2.hist(ass_fit_pars[2], bins=np.arange(min(ass_fit_pars[2]), max(ass_fit_pars[2]) + binwidth, binwidth), edgecolor='black', linewidth=1.2)
    # ax2.set_xlim(0, 1)  # one start
    ax2.set_xlim(0, 3)  # two start
    ax2.set_ylim(0,25)

    # Stacking Energy G1

    binwidth = 1  # kT

    ax3 = fig.add_subplot(2, 2, 2)

    ax3.set_ylabel('count')
    ax3.set_xlabel('G1 (kT)')
    ax3.tick_params(direction='in', top=True, right=True)
    ax3.set_title("Stacking Energy G1")
    ax3.hist(ass_fit_pars[3], bins=np.arange(min(ass_fit_pars[3]), max(ass_fit_pars[3]) + binwidth, binwidth), edgecolor='black', linewidth=1.2)
    ax3.set_xlim(0,35)
    ax3.set_ylim(0,25)

    # G2
    ax4 = fig.add_subplot(2, 2, 4)

    ax4.set_ylabel('count')
    ax4.set_xlabel('G2 (kT)')
    ax4.tick_params(direction='in', top=True, right=True)
    ax4.set_title("Interaction Energy G2")
    ax4.hist(ass_fit_pars[4], bins=np.arange(min(ass_fit_pars[4]), max(ass_fit_pars[4]) + binwidth, binwidth), edgecolor='black', linewidth=1.2)
    ax4.set_xlim(0,35)
    ax4.set_ylim(0,25)

    # # Degeneracy
    # ax5 = fig.add_subplot(2, 3, 6)
    #
    # ax5.set_ylabel('count')
    # ax5.set_xlabel('degeneracy (0..1)')
    # ax5.tick_params(direction='in', top=True, right=True)
    # ax5.set_title("Degeneracy")
    # ax5.hist(ass_fit_pars[5])

    fig.suptitle(p['NRL_str'] + " (n = " + str(len(ass_fit_pars[0])) + ")")

    if show_plot == True:
        plt.show()

    fig.savefig(new_path+p['NRL_str']+"_pars")

    return

def plot_combined_hist(fig, ass_fit_pars, ass_fit_errors, new_path, p, show_plot = True, color='grey', zorder=0, ymax=30):

    ass_fit_pars = np.transpose(np.array(ass_fit_pars))

    # number of stacked nucleosomes in fiber
    ax1 = fig.add_subplot(2, 2, 1)

    ax1.set_ylabel('count')
    ax1.set_xlabel('Stacked Nucleosomes in Fiber')
    ax1.tick_params(direction='in', length=6, width=3, top=True, right=True)
    # ax1.set_title("Stacked Nucleosomes")
    stacked=[]
    for n,x in enumerate(ass_fit_pars[0]):
        stacked.append(ass_fit_pars[0][n]-ass_fit_pars[1][n])
    binwidth = 2 # nucleosomes
    # ax1.hist(stacked, bins=np.arange(min(stacked), max(stacked) + binwidth, binwidth), edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax1.hist(stacked, bins=np.arange(0,35, binwidth), edgecolor='black', linewidth=1.2,
             color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)

    ax1.legend(loc=1, frameon=False)
    ax1.set_xlim(0,40)
    # ax1.set_ylim(-0.1,25)
    ax1.set_ylim(-0.1, ymax)
    # ax1.yaxis.set_ticks(np.arange(0, 30, 5))
    ax1.yaxis.set_ticks(np.arange(0, 35, 5))

    # stiffness
    ax2 = fig.add_subplot(2, 2, 3)

    ax2.set_ylabel('count')
    ax2.set_xlabel('Fiber stiffness (pN/nm)')
    ax2.tick_params(direction='in', length=6, width=3, top=True, right=True)
    # ax2.set_title("Fiber Stiffness")
    # binwidth = 0.15  # pN/nm - 167
    binwidth = 0.05  # pN/nm - 197
    # ax2.hist(ass_fit_pars[2], bins=np.arange(min(ass_fit_pars[2]), max(ass_fit_pars[2]) + binwidth, binwidth), edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax2.hist(ass_fit_pars[2], bins=np.arange(0,3,binwidth),
             edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)

    ax2.legend(loc=1, frameon=False)
    # ax2.set_xlim(-0.03, 1)
    ax2.set_xlim(-0.03, 1.03)
    # ax2.set_ylim(-0.1,25)
    ax2.set_ylim(-0.1,ymax)
    # ax2.yaxis.set_ticks(np.arange(0, 30, 5))
    ax2.yaxis.set_ticks(np.arange(0, 35, 5))

    # Stacking Energy G1

    binwidth = 1.3  # kT

    ax3 = fig.add_subplot(2, 2, 2)

    ax3.set_ylabel('count')
    ax3.set_xlabel('G1 (kT)')
    ax3.tick_params(direction='in', length=6, width=3, top=True, right=True)
    # ax3.set_title("Stacking Energy G1")
    # ax3.hist(ass_fit_pars[3], bins=np.arange(min(ass_fit_pars[3]), max(ass_fit_pars[3]) + binwidth, binwidth), edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax3.hist(ass_fit_pars[3], bins=np.arange(0,30,binwidth),
             edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax3.legend(loc=1, frameon=False)
    # ax3.set_xlim(-0.1,30)
    ax3.set_xlim(-0.1, 35)
    # ax3.set_ylim(-0.1,25)
    ax3.set_ylim(-0.1, ymax)
    # ax3.yaxis.set_ticks(np.arange(0, 30, 5))
    ax3.yaxis.set_ticks(np.arange(0, 35, 5))

    # G2
    ax4 = fig.add_subplot(2, 2, 4)

    ax4.set_ylabel('count')
    ax4.set_xlabel('G2 (kT)')
    ax4.tick_params(direction='in', length=6, width=3, top=True, right=True)
    # ax4.set_title("Interaction Energy G2")
    # ax4.hist(ass_fit_pars[4], bins=np.arange(min(ass_fit_pars[4]), max(ass_fit_pars[4]) + binwidth, binwidth), edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax4.hist(ass_fit_pars[4], bins=np.arange(0,30, binwidth),
             edgecolor='black', linewidth=1.2, color=color, label=p['NRL_str'], zorder=zorder, alpha=0.5)
    ax4.legend(loc=1, frameon=False)
    # ax4.set_xlim(-0.1,30)
    ax4.set_xlim(-0.1, 35)
    # ax4.set_ylim(-0.1,25)
    ax4.set_ylim(-0.1,ymax)
    # ax4.yaxis.set_ticks(np.arange(0, 30, 5))
    ax4.yaxis.set_ticks(np.arange(0, 35, 5))

    return


def fit_circle(X, Y):
    # thanks klaas

    # coordinates of the barycenter
    x_m = np.mean(X)
    y_m = np.mean(Y)

    # calculation of the reduced coordinates
    u = X - x_m
    v = Y - y_m

    # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv = sum(u * v)
    Suu = sum(u ** 2)
    Svv = sum(v ** 2)
    Suuv = sum(u ** 2 * v)
    Suvv = sum(u * v ** 2)
    Suuu = sum(u ** 3)
    Svvv = sum(v ** 3)

    # Solving the linear system
    A = np.array([[Suu, Suv], [Suv, Svv]])
    B = np.array([Suuu + Suvv, Svvv + Suuv]) / 2.0
    uc, vc = np.linalg.solve(A, B)

    xc_1 = x_m + uc
    yc_1 = y_m + vc

    # Calcul des distances au centre (xc_1, yc_1)
    Ri_1 = np.sqrt((X - xc_1) ** 2 + (Y - yc_1) ** 2)
    R_1 = np.mean(Ri_1)
    # residu_1 = sum((Ri_1-R_1)**2)
    return x_m, y_m, R_1


def measurement_pars():
    meas_pars = {}
    meas_pars['date'] = []  # date of measurement
    meas_pars['data'] = []  # measurement
    meas_pars['bead'] = []  # bead
    meas_pars['beads'] = []  # number of beads
    meas_pars['points'] = []  # number of data-points
    meas_pars['points_exp'] = []  # number of expected data-points
    meas_pars['points_frac'] = []  # fraction of data-points
    meas_pars['X0_um'] = []  # global X-position (um)
    meas_pars['Y0_um'] = []  # global Y-position (um)
    meas_pars['Z0_um'] = []  # global Z-position (um)
    meas_pars['dZ_um'] = []  # absolute Z-extension after drift correction(um)
    meas_pars['amp'] = []  # mean amplitude
    meas_pars['drift_nm_s'] = []  # drift (nm/s)
    meas_pars['L_bp'] = []  # number of base pairs
    meas_pars['P_nm'] = []  # persistence length (nm)
    meas_pars['S_pN'] = []  # stretch modulus (pN)
    meas_pars['NRL'] = []  # nucleosome repeat length
    meas_pars['repeats'] = []  # number of repeats
    meas_pars['N_nuc'] = []  # number of steps
    meas_pars['N_tet'] = []  # number of tetrasomes
    meas_pars['N_stack'] = []  # number of stacked nucleosomes
    meas_pars['k_pN_nm'] = []  # fiber stiffness (pN/nm)
    meas_pars['k_pN_nm_se'] = []  # standard error fiber stiffness (pN/nm)
    meas_pars['G1_kT'] = []  # G1 (kT)
    meas_pars['G1_kT_se'] = []  # standard error G1 (kT)
    meas_pars['G2_kT'] = []  # G2 (kT)
    meas_pars['G2_kT_se'] = []  # standard error G2 (kT)
    meas_pars['G3_kT'] = []  # G3 (kT)
    meas_pars['L_fold'] = []  # folded length (nm/nuc) - standard: 1.5 nm/nuc
    meas_pars['L_unwrap'] = []  # unwrapped length (basepairs/nucleosome) - standard: 56 bp/nuc
    meas_pars['L_ext'] = []  # extended length (nm/nucleosome) - standard: 5 nm/nuc
    meas_pars['degeneracy'] = []  # degeneracy (0 = two-start, 1 = one-start helix)
    meas_pars['f_rupt_pN'] = []  # rupture force (pN)
    meas_pars['radius_um'] = []  # radius of circle (um)
    return meas_pars