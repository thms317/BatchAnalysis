import numpy as np
import pandas as pd
import os
import glob
import functions as func
import shutil

folder = "C:\\Users\\brouw\\Desktop\\167 x 30\\"
recalc_folder = "C:\\Users\\brouw\\Desktop\\167 x 30\\recalculated\\"


# recalculate fitfiles + save

fitfiles = []
os.chdir(folder)
for file in glob.glob("*.fit"):
    fitfiles.append(file)

for file in fitfiles:

    # read DataFrame
    df = pd.read_csv(folder+file, sep="\t")

    force = np.array(df['F (pN)'])
    time = np.array(df['t (s)'])

    A = 59.5
    C = 0.01
    L = 1.4

    magnet = - np.log(force/A - C/A) * L

    recalc_force = func.calc_force(magnet)

    df['F (pN)'] = recalc_force

    df.to_csv(recalc_folder+file[:-4]+"_rec.fit", sep='\t')


# rename logfiles + save

logfiles = []
os.chdir(folder)
for file in glob.glob("*.log"):
    logfiles.append(file)

for file in logfiles:

    old_name = folder+file

    base, extension = os.path.splitext(file)

    new_name = recalc_folder + base + "_rec.log"

    shutil.copy(old_name, new_name)

