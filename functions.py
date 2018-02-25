import numpy as np
from scipy.signal import argrelextrema
from scipy import stats
import mpmath

# worm-like chain
def WLC(f, p, L, S, x0):
    kBT = 4.114  # (pN nm) - Boltzmann factor
    return (L * (1 - 0.5 * (np.sqrt(kBT / (f * p))) + f / S)) / 1000 + x0

# Langevin equation (needed for FJC)
def Langevin(x):
    return mpmath.coth(x) - 1 / x


# Freely Jointed Chain
def FJC(f, b):
    kBT = 4.114  # (pN nm) - Boltzmann factor
    return L * Langevin((f * b) / kBT) + L * f / S


# # FJC by John
# def FJC_john(f, b=None, k_pN_nm=0.1, L_nm=20, S_pN=1e3):
#     kBT = 4.114  # (pN nm) - Boltzmann factor
#     if b == None:
#         b = 3 * kBT / (k_pN_nm * L_nm)
#     x = f * b / kBT
#     z = L_nm * (sympy.coth(x) - 1 / x)
#     z = np.asarray(z)
#     z += L_nm * f / S_pN
#     return z

# L_app
def NewP(L, p, alpha, i):
    C = 8 * (1 - np.cos(alpha / 4))
    ro = i / L
    return p / ((1 + p * i * C * ro) ** 2)


# degrees to radians
def degrad(deg):
    return deg * (np.pi / 180)


# calculate force
def calc_force(i):
    A = 85  # for 2.8 um beads (pN)
    l1 = 1.4  # decay length 1 (mm)
    l2 = 0.8  # decay length 2 (mm)
    f0 = 0.01  # force-offset (pN)
    return A * (0.7 * np.exp(-i / l1) + 0.3 * np.exp(-i / l2)) + f0


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


# function to reject outliers, replace with NaN
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


# get numbers from string
def get_num(x):
    return float(''.join(ele for ele in x if ele.isdigit() or ele == '.'))


# find peaks
def peak_finder_old(x, y):
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


def peak_finder(y):  # Finds y peaks at position x in xy graph
    y = np.array(y)
    # x=list(x)

    # mirror the data in the Y-axis (to find potential peak at x=0)
    # x_x = list(reversed(np.negative(np.array(x[1:])))) + x
    Yy = np.append(y[:-1], y[::-1])
    yYy = np.append(y[::-1][:-1], Yy)

    from scipy.signal import argrelextrema
    maxInd = argrelextrema(np.array(yYy), np.greater)
    r = np.array(yYy)[maxInd]
    a = maxInd[0]

    # discard all peaks for negative dimers
    peaks_index = []
    peaks_height = []
    for n, i in enumerate(a):
        i = 1 + i - len(y)
        if i >= 0 and i <= len(y):
            peaks_height.append(r[n])
            peaks_index.append(i)

    return peaks_index, peaks_height


# fitting two lines
def two_lines(x, a, b, c, d):
    one = a * x + b
    two = c * x + d
    return np.maximum(one, two)


# calculate drift using stuck bead
def calc_drift_stuck(data_lines, headers, time, beads):
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


# calculate drift using self
def calc_drift_self(data_lines, headers, time, bead):
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