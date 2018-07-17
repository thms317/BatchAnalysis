import pandas as pd
import numpy as np
import os
import itertools as it

folder = "S:\\Brouwer\\Chromatin Force Spectroscopy\\Cummulative\\assembled\\"

titles, k, G1, G2 = [], [], [], []

os.chdir(folder)
for n, pars_file in enumerate(os.listdir(os.getcwd())):
    title = pars_file[:6]
    titles.append(title)

    df = pd.read_csv(folder+pars_file, sep='\t')

    k.append(np.array(df['k_pN_nm']))
    G1.append(np.array(df['G1_kT']))
    G2.append(np.array(df['G2_kT']))

k_pars = list(it.zip_longest(*k, fillvalue = ''))
G1_pars = list(it.zip_longest(*G1, fillvalue = ''))
G2_pars = list(it.zip_longest(*G2, fillvalue = ''))

saveloc = "C:\\Users\\brouw\\Desktop\\"

df1 = pd.DataFrame(data=k_pars)
df1.columns = titles
df1.to_csv(saveloc + "k_pars.dat", sep='\t')

df2 = pd.DataFrame(data=G1_pars)
df2.columns = titles
df2.to_csv(saveloc + "G1_pars.dat", sep='\t')

df3 = pd.DataFrame(data=G2_pars)
df3.columns = titles
df3.to_csv(saveloc + "G2_pars.dat", sep='\t')

