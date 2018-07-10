import numpy as np
from scipy import stats

kT = 4.114  # pN nm
L0 = 0.34  # nm / bp


# worm-like chain
def WLC(f, L_bp=1e3, P_nm=50, S_pN=1e3, z0_nm=0, HMf_cWLC=0):
    z = 1 - 0.5 * np.sqrt(kT / (f * P_nm)) + f / S_pN
    z *= L_bp * 0.34
    z += z0_nm
    # w = f dz= f*z - z df
    w = f * z - L_bp * L0 * (f - np.sqrt((f * kT) / P_nm) + f ** 2 / (2 * S_pN))
    return z, w / kT


def WLC_fit(f, p, L, S, z0):
    kT = 4.114  # (pN nm) - Boltzmann factor
    return (L * (1 - 0.5 * (np.sqrt(kT / (f * p))) + f / S)) / 1000 + z0


# freely jointed chain
def FJC(f, b=None, k_pN_nm=0.01, L_nm=5, z0_nm=0, HMf_cFJC=0):
    if b == None:
        b = 3 * kT / (k_pN_nm * L_nm)
    x = f * b / kT
    # coth(x)= (exp(x) + exp(-x)) / (exp(x) - exp(-x)) --> see Wikipedia
    exp_x = np.exp(x)
    z = (exp_x + (1 / exp_x)) / (exp_x - (1 / exp_x)) - 1 / x
    z *= L_nm
    z += z0_nm
    # w = f dz= f*z - z df
    # z = coth(x) - 1/x
    z_df = L_nm * (kT / b) * (np.log(np.sinh(x)) - np.log(x))
    # z_df = L_nm * (kT / b) * (np.log(np.sinh(x)) - np.log(x)) # 4 pi missing
    w = f * z - z_df
    return z, w / kT


# L_app - function modified
def kCalcP(L_bp, P_nm=50, N=0, alpha=0):
    #   Persistence length correction according to Kulic and Schiessel
    C = 8 * (1 - np.cos(alpha * (np.pi / 180) / 4.0))
    ro = N / (L_bp * L0)
    return P_nm / ((1 + P_nm * N * C * ro) ** 2)


# degrees to radians
def degrad(deg):
    return deg * (np.pi / 180)


# radians to degrees
def raddeg(rad):
    return rad * (180 / np.pi)


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

def get_int(string):
    return int(''.join(x for x in string if x.isdigit()))

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
def drift_stuck(data_lines, headers, time, beads):
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
def drift_self(Z, time):
    # correcting the drift using a stuck bead
    Z_start = np.percentile(Z[:100] * 1000, 1)
    Z_end = np.percentile(Z[-100:] * 1000, 1)
    t_start = np.percentile(time[:100], 1)
    t_end = np.percentile(time[-100:], 1)
    drift_slope, intercept, r_value, p_value, std_err = stats.linregress([t_start, t_end], [Z_start, Z_end])
    return drift_slope


# calculate ruptures
def rupture(time, amplitude, mask=False):
    # calculating the first derivative of amplitude
    dx = np.diff(time)
    dy = np.diff(amplitude)
    diff_amp = np.append([0], np.divide(dy, dx))

    diff_amp -= np.mean(diff_amp)
    std_diff_amp = abs(np.std(diff_amp))
    peak = max(abs(diff_amp))


    # classify as 'tether break' if absolute diff(amplitude) maximum exceeds n times the std
    n = 10

    if mask:
        mask = np.ones(len(amplitude))
        if peak > n * std_diff_amp:
            peak_index = int(np.where(abs(diff_amp) == peak)[0])
            average_before = np.mean(amplitude[:peak_index])
            average_after = np.mean(amplitude[peak_index:])
            mean_diff = abs(average_after - average_before)
            if mean_diff > 2 * abs(np.std(amplitude)):
                mask_on = np.ones(peak_index)
                fuck_it_mask_off = np.zeros(len(amplitude) - peak_index)
                mask = np.concatenate((mask_on, fuck_it_mask_off))
                return True, peak_index, mask
            else:
                return False, peak_index, mask
        else:
            peak_index = None
            return False, peak_index, mask


    else:
        if peak > n * std_diff_amp:
            peak_index = int(np.where(abs(diff_amp) == peak)[0])
            average_before = np.mean(amplitude[:peak_index])
            average_after = np.mean(amplitude[peak_index:])
            mean_diff = abs(average_after - average_before)
            if mean_diff > 2 * abs(np.std(amplitude)):
                return True, peak_index
            else:
                return False, peak_index
        else:
            peak_index = None
            return False, peak_index