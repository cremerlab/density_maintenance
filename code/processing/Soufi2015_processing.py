#%%
import numpy as np 
import pandas as pd 


# Load the data and dropped unnecessary columns
data = pd.read_csv('../../data/source/Soufi2015/soufi2015_longform_annotated.csv')
data.drop(columns=['Unnamed: 0', 'reported_tot_per_cell', 'corrected_volume', 'dataset'],
                inplace=True)

# Compute the mass fractions
dfs = []
for g, d in data.groupby(['growth_rate_hr']):
    tot = d['reported_fg_per_cell'].sum()
    d['mass_frac'] = d['reported_fg_per_cell'].values / tot
    dfs.append(d)

# Concatenate the dataset and drop an unnecessary column
data = pd.concat(dfs, sort=False)
data.drop(columns=['reported_fg_per_cell'], inplace=True)

# Save
data.to_csv('../../data/source/Soufi2015/Soufi2015_processed.csv', index=False)