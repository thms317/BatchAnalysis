import pandas as pd
import numpy as np
import functions as func
import matplotlib.pyplot as plt

measurements = []
measurements.append(["168x16", "180612", "data_028", "34"])
# measurements.append(["169x16", "180621", "data_036", "54"])  # 1/2
measurements.append(["169x16", "180621", "data_062", "36"])  # 2/2
measurements.append(["171x16", "180629", "data_053", "22"])
measurements.append(["173x16", "180706", "data_010", "0"])  # 1/2
# measurements.append(["173x16", "180705", "data_033", "21"])  # 2/2
# measurements.append(["175x16", "180712", "data_033", "20"])  # 1/3
# measurements.append(["175x16", "180712", "data_021", "48"])  # 2/3
measurements.append(["175x16", "180712", "data_015", "9"])  # 3/3
measurements.append(["177x16", "180712", "data_052", "21"])

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Selected measurements\\"
pars_file = "S:\\Brouwer\\Chromatin Force Spectroscopy\\all_pars.dat"

plt.figure(figsize=(12,6))

for n, measurement in enumerate(measurements):
    file = measurement[1] + "_" + measurement[2] + "_" + measurement[3] + ".fit"
    df = pd.read_csv(folder + file, sep="\t")

    # columns
    time = np.array(df['t (s)'])
    force = np.array(df['F (pN)'])
    z = np.array(df['z (um)'])
    z_fit = np.array(df['z fit (um)'])

    # get rupture force
    df_pars = pd.read_csv(pars_file, sep="\t")
    f_rupt = df_pars.loc[
        (df_pars['date'] == int(measurement[1])) & (df_pars['data'] == func.get_int((measurement[2]))) & (
                    df_pars['bead'] == float(measurement[3])), "f_rupt_pN"]

    # split in pull & release
    dx = np.diff(time)
    dy = np.diff(force)
    diff_force = np.append([0], np.divide(dy, dx))  # add a zero as first element

    factor = max(diff_force / 1000)

    trajectory = "T1"

    if measurement[0] == "175x16" or measurement[0] == "177x16" :
        trajectory = "T2"


    if trajectory == "T1":
        f_pull = force[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        z_pull = z[np.where((diff_force > factor) & (time > 50) & (time < 80))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 50) & (time < 80))]
    if trajectory == "T2":
        f_pull = force[np.where((diff_force > factor) & (time > 90) & (time < 125))]
        z_pull = z[np.where((diff_force > factor) & (time > 90) & (time < 125))]
        z_fit_pull = z_fit[np.where((diff_force > factor) & (time > 90) & (time < 125))]

    offset = 0.3
    z_pull += (n*offset)
    z_fit_pull += (n*offset)

    plt.scatter(z_pull,f_pull,label=measurement[0], s = 5)
    # plt.plot(z_fit_pull, f_pull, color ="black", linewidth = 2)

plt.xlim(0.5,3.25)

plt.legend()
plt.title("Force Spectoscopy various NRLs")

plt.savefig("C:\\Users\\brouw\\Desktop\\combined.png",dpi=300)
plt.xlim(0.5,3)
plt.ylim(0,10)
plt.savefig("C:\\Users\\brouw\\Desktop\\combined_zoom.png",dpi=300)

plt.close()

