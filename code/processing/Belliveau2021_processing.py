#%%
import numpy as np 
import pandas as pd 

data = pd.read_csv('../../data/source/Belliveau2021/compiled_absolute_measurements.csv')

# Compute the mass fractions
dfs = []
for g, d in data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    tot = d['fg_per_cell'].sum()
    frac = d['fg_per_cell'].values / tot
    d['mass_frac'] = frac
    dfs.append(d)
data = pd.concat(dfs, sort=False)

# Drop unnecessary columns
data.drop(columns=['gene_product', 'tot_per_cell', 'fg_per_cell', 'dataset'],
            inplace=True)
data.to_csv('../../data/source/Belliveau2021/Belliveau2021_processed.csv', index=False)
# %%
