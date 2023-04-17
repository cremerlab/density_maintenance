# %%
import numpy as np
import pandas as pd
import tqdm
babu = pd.read_csv(
    '../../data/literature/Babu2018/Babu2018_minimal_classification.csv')
files = ['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']
dfs = [pd.read_csv(
    f'../../data/literature/{f}/{f}_processed.csv') for f in files]
data = pd.concat(dfs, sort=False)
dfs = []
for g, d in data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    if np.round(d['mass_frac'].sum(), decimals=3) == 1.0:
        dfs.append(d)
data = pd.concat(dfs, sort=False)

data.drop(columns=['cog_desc', 'cog_category', 'gene_product', 'annotation', 'dataset'],
          inplace=True)

# babu = babu[babu['localization'] != 'EC']
data['envelope'] = False
data['periplasm'] = False
data['membrane'] = False
data['ribosomal'] = False
data['metabolism'] = False
for g, _ in tqdm.tqdm(babu.groupby(['name', 'b_number', 'localization'])):
    _dict = {'periplasm': False, 'envelope': True, 'membrane': False}
    if g[-1] == 'PE':
        _dict['periplasm'] = True
    else:
        _dict['periplasm'] = False
    if g[-1] in ['IM', 'OM', 'LPO', 'LPI', 'MR']:
        _dict['membrane'] = True
    else:
        _dict['membrane'] = False

    if (len(data[data['gene_name'] == g[0]]) == 0) & (len(data[data['b_number'] == g[1]]) == 0):
        continue

    for k, v in _dict.items():
        data.loc[data['gene_name'] == g[0], k] = v
        data.loc[data['b_number'] == g[1], k] = v

data.loc[data['cog_letter'].isin(
    ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q']), 'metabolism'] = True
data.loc[data['cog_letter'].isin(['J']), 'ribosomal'] = True
data = data[data['dataset_name'] != 'Valgepea et al. 2013']
data.to_csv('../../data/literature/compiled_mass_fractions.csv', index=False)

# %%

# %%

# %%

# %%

# %%
