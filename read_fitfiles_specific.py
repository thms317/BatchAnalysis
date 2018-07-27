import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functions as func

folder = "C:\\Users\\tbrouwer\\Desktop\\Parameters vs Number of Repeats\\167x15 + 197x15\\"
save_folder = "C:\\Users\\tbrouwer\\Desktop\\Parameters vs Number of Repeats\\"

measurements = []
# measurements.append(["160121_data_010_53.fit", "197x15"])
measurements.append(["160121_data_010_60.fit", "197x15"])  # best
measurements.append(["15x167 FC1_data_006_40.fit", "167x15"])

fig = plt.figure(figsize=(30, 10))
plt.rcParams.update({'font.size': 30})  # legend size
plt.rc('axes', linewidth=3)
plt.rc('xtick', labelsize=30)
plt.rc('ytick', labelsize=30)

for measurement in measurements:

    # open data
    file_all = folder + measurement[0]

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

    # calculating the first derivative of force
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element
    factor = max(diff_force / 1000)

    if measurement[1] == "167x15":
        f_pull = force[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        f_release = force[np.where((diff_force < factor) & (time > 140))]
        z_pull = z[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        z_release = z[np.where((diff_force < factor) & (time > 140))]
        time_pull = time[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        time_release = time[np.where((diff_force < factor) & (time > 140))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        T1 = T1[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        T2 = T2[np.where((diff_force > factor) & (time > 115) & (time < 155))]
        T3 = T3[np.where((diff_force > factor) & (time > 115) & (time < 155))]
    else:
        f_pull = force[np.where((diff_force > factor))]
        f_release = force[np.where((diff_force < factor))]
        z_pull = z[np.where((diff_force > factor))]
        z_release = z[np.where((diff_force < factor))]
        time_pull = time[np.where(diff_force > factor)]
        time_release = time[np.where((diff_force < factor))]
        z_fit_pull = z_fit[np.where((diff_force > factor))]
        T1 = T1[np.where((diff_force > factor))]
        T2 = T2[np.where((diff_force > factor))]
        T3 = T3[np.where((diff_force > factor))]

    transitions = np.stack((T1, T2, T3))  # all transitions in a 3D array

    if measurement[1] == "197x15":

        # wlc
        f_wlc = np.linspace(0.1,35, 1000)
        z_wlc, _ = func.WLC(f_wlc, L_bp=4985, S_pN=900)

        ax1 = fig.add_subplot(1, 2, 2)
        ax1.scatter(z_pull, f_pull, label=measurement[1], color='blue', s=100, zorder=25, facecolors='none', alpha=0.5)
        for t in range(len(np.transpose(transitions[0]))):
            plt.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3
        # ax1.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none', alpha=0.5)
        ax1.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)
        ax1.plot(z_wlc/1000, f_wlc, color='black', linewidth=3, zorder=10)

        # ax1.legend(loc=2, frameon=False)
        ax1.set_xlim(-0.1, 1.75)
        ax1.set_ylim(-1, 35)
        ax1.tick_params(direction='in', length=6, width=3, top=True, right=True)
        ax1.xaxis.set_ticks(np.arange(0, 2, 0.5))
        ax1.yaxis.set_ticks(np.arange(0, 40, 10))

        ax1.set_ylabel('F (pN)')
        ax1.set_xlabel('z ($\mu$m)')
        ax1.set_title(measurement[1])

    elif measurement[1] == "167x15":

        # wlc
        f_wlc = np.linspace(0.1,35, 1000)
        z_wlc, _ = func.WLC(f_wlc, L_bp=4500, S_pN=900)

        ax2 = fig.add_subplot(1, 2, 1)
        ax2.scatter(z_pull, f_pull, label=measurement[1], color='red', s=100, zorder=25, facecolors='none', alpha=0.5)
        for t in range(len(np.transpose(transitions[0]))):
            plt.plot(np.transpose(transitions[2])[t], f_pull, '--', color='lightgrey')  # transition 3
        # ax2.scatter(z_release, f_release, color='lightgrey', s=30, zorder=15, facecolors='none', alpha=0.5)
        ax2.plot(z_fit_pull, f_pull, color='black', linewidth=3, zorder=1000)
        ax2.plot(z_wlc / 1000, f_wlc, color='black', linewidth=3, zorder=10)

        # ax2.legend(loc=2, frameon=False)
        ax2.set_xlim(-0.1, 1.75)
        ax2.set_ylim(-1, 35)
        ax2.tick_params(direction='in', length=6, width=3, top=True, right=True)
        ax2.xaxis.set_ticks(np.arange(0, 2, 0.5))
        ax2.yaxis.set_ticks(np.arange(0, 40, 10))

        ax2.set_ylabel('F (pN)')
        ax2.set_xlabel('z ($\mu$m)')
        ax2.set_title(measurement[1])

plt.savefig(save_folder+"167x15 + 197x15.png", dpi=600)
# plt.show()
