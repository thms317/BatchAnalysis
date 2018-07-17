import pandas as pd
import os

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\assembled\\"
saveloc = "S:\\Brouwer\\Chromatin Force Spectroscopy\\"

os.chdir(folder)
for n, file in enumerate(os.listdir(os.getcwd())):
    if n == 0:
        df = pd.read_csv(folder + file, sep="\t")
    else:
        df2 = pd.read_csv(folder + file, sep="\t")
        df = df.append(df2)

df.to_csv(saveloc + "all_pars.dat", sep='\t')