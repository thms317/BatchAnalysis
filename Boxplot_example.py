import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

datafile = "PlantGrowth.csv"
folder = "C:\\Users\\tbrouwer\\Downloads\\"
data = pd.read_csv(folder + datafile)

# Create a boxplot
data.boxplot('weight', by='group', figsize=(12, 8))

mod = ols('weight ~ group', data = data).fit()
aov_table = sm.stats.anova_lm(mod, type = 2)
print(aov_table)

ctrl = data['weight'][data.group == 'ctrl']

grps = pd.unique(data.group.values)
d_data = {grp: data['weight'][data.group == grp] for grp in grps}

k = len(pd.unique(data.group))  # number of conditions
N = len(data.values)  # conditions times participants
n = data.groupby('group').size()[0]  # Participants in each condition

plt.show()
