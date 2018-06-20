import pandas as pd
import numpy as np
import functions as func
import matplotlib.pyplot as plt

good = "C:\\Users\\brouw\\Desktop\\13_good_21_bad\\data_013.dat"
bad = "C:\\Users\\brouw\\Desktop\\13_good_21_bad\\data_021.dat"

# read DataFrame
# df = pd.read_csv(good, sep="\t")
df = pd.read_csv(bad, sep="\t")

time = df['Time (s)']
magnet = df['Stepper shift (mm)']

# for n in range(len(df.columns)):
#     column = np.array(df.ix[:,n])
#     bool = np.isnan(column)
#     if True in bool:
#         print("Something is amiss")

force = func.calc_force(magnet)

# print(np.isnan(force).any())

# calculating the first derivative of magnet
dx = np.diff(time)
dy = np.diff(magnet)
# diff_magnet = np.append([0], np.divide(dy, dx))  # add a zero as first element

plt.scatter(np.arange(0,len(dx)),dx)
plt.show()