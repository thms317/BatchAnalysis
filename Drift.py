def calc_drift(data_lines,headers,time,beads):

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

    # actually correcting drift
    Z_drift = []
    for n, t in enumerate(time):
        Z_drift.append(Z[n] - (slope / 1000) * t)
    Z = Z_drift