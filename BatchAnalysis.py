# import relevant stuff
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.stats import norm
from scipy import stats
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.signal import argrelextrema
import os

# definitions

# L_app
def NewP(L, p, alpha, i):
    C = 8 * (1 - np.cos(alpha / 4))
    ro = i / L
    return p / ((1 + p * i * C * ro) ** 2)

# degrees to radians
def degrad(deg):
    return deg * (np.pi / 180)

# calculate force from magnet position
def calc_force(i):
    A = 85  # for 2.8 um beads (pN)
    l1 = 1.4  # decay length 1 (mm)
    l2 = 0.8  # decay length 2 (mm)
    f0 = 0.01  # force-offset (pN)
    return A * (0.7 * np.exp(-i / l1) + 0.3 * np.exp(-i / l2)) + f0

# worm-like chain
def WLC(f, p, L, S, x0):
    return (L * (1 - 0.5 * (np.sqrt(kBT / (f * p))) + f / S)) / 1000 + x0

# median filter
def medfilt(x, k):  # Apply a length-k median filter to a 1D array x. Boundaries are extended by repeating endpoints.
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros((len(x), k), dtype=x.dtype)
    y[:, k2] = x
    for i in range(k2):
        j = k2 - i
        y[j:, i] = x[:-j]
        y[:j, i] = x[0]
        y[:-j, -(i + 1)] = x[j:]
        y[-j:, -(i + 1)] = x[-1]
    return np.median(y, axis=1)

# reject outliers, replace with NaN
def reject_outliers(data):
    data_filtered = []
    norm_data = []
    norm_data = abs(data - np.mean(data))
    for n, i in enumerate(norm_data):
        if i > 2 * np.std(data):
            i = np.nan
        else:
            i = data[n]
        data_filtered.append(i)
    return data_filtered

# isolate numbers from string/path
def get_num(x):
    return float(''.join(ele for ele in x if ele.isdigit() or ele == '.'))

# find peaks in data
def peak_finder(x, y):
    x = list(x)
    y = list(y)

    # mirror the data in the Y-axis (to find potential peak at x=0)

    x_x = list(reversed(np.negative(np.array(x[1:])))) + x
    y_y = list(reversed(y[1:])) + y

    maxInd = argrelextrema(np.array(y_y), np.greater)
    r = np.array(y_y)[maxInd]
    a = maxInd[0]

    # discard all peaks for negative dimers

    peaks_index = []
    peaks_height = []
    for n, i in enumerate(a):
        i = 1 + i - len(y)
        if i >= 0:
            peaks_height.append(r[n])
            peaks_index.append(x[i])

    return peaks_index, peaks_height

# calculate drift
def calc_drift_stuck(data_lines,headers,time,beads):

    # find a stuck bead by RMS analysis of single bead in Z-direction
    z_rms = []
    for b in range(0, beads):
        z_temp = []
        for x in data_lines:
            z_temp.append(float(x.split()[headers.index('Z' + str(b) + ' (um)')]))
        z_temp = np.array(z_temp)
        z_temp -= np.mean(z_temp)
        z_rms.append(np.sqrt(np.mean(z_temp ** 2)))
    stuck_index = int(z_rms.index(min(z_rms)))

    z_stuck = []
    for x in data_lines:
        z_stuck.append(float(x.split()[headers.index('Z' + str(stuck_index) + ' (um)')]))

    # correcting the drift using a stuck bead
    driftz = []
    driftt = []
    minmax = []
    for n, z in enumerate(z_stuck):
        driftt.append(time[n])
        driftz.append(z * 1000)
    minmax.append(np.percentile(driftz[:100], 1))
    minmax.append(np.percentile(driftz[-100:], 1))
    minmax.append(np.percentile(driftt[:100], 1))
    minmax.append(np.percentile(driftt[-100:], 1))
    slope, intercept, r_value, p_value, std_err = stats.linregress([minmax[2], minmax[3]], [minmax[0], minmax[1]])

    return slope

def calc_drift_self(data_lines,headers,time,bead):

    Z = []
    for x in data_lines:
        Z.append(float(x.split()[headers.index('Z' + str(bead) + ' (um)')]))

    # correcting the drift using a stuck bead
    driftz = []
    driftt = []
    minmax = []
    for n, z in enumerate(Z):
        driftt.append(time[n])
        driftz.append(z * 1000)
    minmax.append(np.percentile(driftz[:100], 1))
    minmax.append(np.percentile(driftz[-100:], 1))
    minmax.append(np.percentile(driftt[:100], 1))
    minmax.append(np.percentile(driftt[-100:], 1))
    slope, intercept, r_value, p_value, std_err = stats.linregress([minmax[2], minmax[3]], [minmax[0], minmax[1]])

    return slope


# constants

kBT = 4.114  # (pN nm) - Boltzmann factor
Lc = 500  # contour length (bp)
L = Lc * 0.34  # contour length (nm)
p = 50  # persistence length (nm)
S = 1000  # stretch modulus (pN)
x0 = 0  # offset (nm)


# paths
data_path = "C:\\Users\\brouw\\Desktop\\"
subfolder = "DNA mono 601 - 5-10"
analysis_path = data_path + subfolder+"\\"
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)

# headers
sub = []
measurements = []

# measurements
# DNA
measurements.append(['180201', '007'])
measurements.append(['171206', '003'])
measurements.append(['171206', '005'])
measurements.append(['171206', '006'])
measurements.append(['171206', '007'])
measurements.append(['171206', '008'])
measurements.append(['171206', '009'])
measurements.append(['171206', '010'])
measurements.append(['171206', '011'])
measurements.append(['171206', '012'])

measurements.append(['171206', '022'])
measurements.append(['171206', '023'])
measurements.append(['171206', '024'])
measurements.append(['171206', '012'])
measurements.append(['171206', '025'])
measurements.append(['171206', '026'])
measurements.append(['171206', '027'])
measurements.append(['171206', '028'])
measurements.append(['171206', '029'])

# chromatin
#measurements.append(['171206', '013'])
#measurements.append(['171206', '014'])
#measurements.append(['171206', '016'])
#measurements.append(['171206', '017'])
#measurements.append(['171206', '018'])
#measurements.append(['171206', '019'])
#measurements.append(['171206', '020'])
#measurements.append(['171206', '021'])


for n, i in enumerate(measurements):

    # open data
    sub = measurements[n]
    file_location = data_path + str(sub[0]) + "\\"
    file_name = "data_" + str(sub[1])
    file_extension = ".dat"
    file_all = file_location + file_name + file_extension

    # title - extract date of the measurements
    title = int(get_num(file_location))

    # headers
    f = open(file_all, 'r')
    headers_lines = f.readlines()[0]
    f.close()

    headers = headers_lines.split('\t')

    # data
    f = open(file_all, 'r')
    data_lines = f.readlines()[1:]
    f.close()

    # number of beads
    file_name = "data_" + str(sub[1])
    file_extension = ".log"
    file_all = file_location + file_name + file_extension

    # find "number of beads" column
    f = open(file_all, 'r')
    beads = f.readlines()[9]
    f.close()
    beads = int(get_num(beads))

    # calculate force from magnet position
    magnet = []
    time = []
    force = []
    for n, x in enumerate(data_lines):
        magnet.append(float(x.split()[headers.index('Stepper shift (mm)')]))
        time.append(float(x.split()[headers.index('Time (s)')]))
    magnet = np.array(magnet)
    for i in magnet:
        force.append(calc_force(i))

    # calculating the first derivative of magnet to discriminate in pull/release-curve
    dx = np.diff(time)
    dy = np.diff(magnet)
    diff_magnet = np.array(np.divide(dy, dx))

#    # calculate drift for dataset
#    slope = calc_drift_stuck(data_lines, headers, time, beads)

    for bead in range(0, beads):
        
        # load the data
        Z = []
        for x in data_lines:
            Z.append(float(x.split()[headers.index('Z' + str(bead) + ' (um)')]))
        
        # calculate drift for the individual bead
        slope = calc_drift_self(data_lines, headers, time, bead)
        
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
        time_pull = []
        time_release = []
        trigger = []  # from what data point does the pulling trace start

        # if the differential of the magnet is positive -> pull, else -> release ('factor' since 0 does not work)
        for n, i in enumerate(diff_magnet):
            factor = max(diff_magnet / 1000)
            if i < -factor:
                trigger.append(n)  # pull the trigger
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
            wlc.append(WLC(f, p, L, S, x0))

        # select data
        select_f = []
        select_z = []
        for n, f in enumerate(f_pull):
            if 5 < f < 10:
                select_f.append(f)
                select_z.append(Z[n + min(trigger)])

        # initial guesses
        x_init = 1

        # fit the WLC in fashion (x,y) - only fit offset, fix everything else
        popt, pcov = curve_fit(lambda f, x0: WLC(f, p, L, S, x0), select_f, select_z, p0=(x_init))
        std = np.sqrt(np.diag(pcov))  # returns the standard deviation

        x_fit = popt[0]

        z_pull -= x_fit  # subtract fitted offset from data
        z_release -= x_fit  # subtract fitted offset from data
        select_z -= x_fit

        # seperate plots
        plt.title(str(title) + " / " + str(file_name) + " / bead " + str(bead))
        plt.ylabel("Force (pN)")
        plt.xlabel("Extension ($\mu$m)")
        plt.xlim(0, 0.5)
        plt.plot(wlc, f_pull, color='black', zorder=100)
        plt.scatter(z_pull, f_pull, color='darkgreen')
        plt.savefig(analysis_path + str(title) + "_" + str(file_name) + "_bead" + str(bead) + '.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        plt.title(str(title) + " / " + str(file_name) + " / bead " + str(bead))
        plt.xlabel("Time (s)")
        plt.ylabel("Extension ($\mu$m)")
        plt.scatter(time, Z, color='darkgreen')
        plt.savefig(analysis_path + str(title) + "_" + str(file_name) + "_bead" + str(bead) + '_time.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # subplots
        plt.title(str(title) + " / " + str(file_name) + " / bead " + str(bead))

        plt.subplot(2, 1, 1)
        plt.plot(wlc, f_pull, color='black', zorder=100)
        plt.scatter(z_pull, f_pull, color='darkgreen')
        plt.ylabel("Force (pN)")
        plt.xlabel("Extension ($\mu$m)")
        plt.xlim(0, 0.5)
        
        plt.subplot(2, 1, 2)
        plt.scatter(time, Z, color='darkgreen')
        plt.xlabel("Time (s)")
        plt.ylabel("Extension ($\mu$m)")
        
        plt.savefig(analysis_path + str(title) + "_" + str(file_name) + "_bead" + str(bead) + '_subplot.png', dpi=300, bbox_inches='tight')
        
        plt.show()

