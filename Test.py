import json
import pandas as pd

dir = "C:\\Users\\tbrouwer\\Desktop\\testing\\"
json_file = "172x12_human_assembled_pars.json"
df_file = "dataframe.dat"

with open(dir + json_file, 'r') as read_file:
    list_of_dict = json.loads(read_file.read())

values = []
for n, dict in enumerate(list_of_dict):
    if n == 0:
        keys = list(dict.keys())
    values.append(list(dict.values()))

df = pd.DataFrame(data=values)
df.columns = keys

df.to_csv(dir + df_file, sep='\t')