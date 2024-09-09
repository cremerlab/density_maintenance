# %%
import numpy as np
import pandas as pd
import tqdm
babu = pd.read_csv(
    '../../../data/literature/Babu2018/Babu2018_minimal_classification.csv')
files = ['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']  
dfs = [pd.read_csv(
    f'../../../data/literature/{f}/{f}_processed.csv') for f in files]
data = pd.concat(dfs, sort=False)

# Exclude valgepea et al. which is not well designed to measure membrane proteins
data = data[data['dataset_name'] != 'Valgepea et al. 2013']


#%%
# %%
dfs = []
for g, d in data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    if np.round(d['mass_frac'].sum(), decimals=3) <= 1.0:
        dfs.append(d)
data = pd.concat(dfs, sort=False)
data = data[data['dataset_name'] != 'Valgepea et al. 2013']
# data = data[data['dataset_name'] != 'Mori et al. 2021']
data.drop(columns=['cog_desc', 'cog_category', 'gene_product', 'annotation', 'dataset'],
          inplace=True)
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']
# %%
# babu = babu[babu['localization'] != 'EC']
data['envelope'] = False
data['periplasm'] = False
data['membrane'] = False
data['ribosomal'] = False
data['metabolism'] = False
data['inner membrane'] = False
data['outer membrane'] = False
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
    if g[-1] in ['IM', 'LPI']:
        _dict['inner membrane'] = True
    else:
        _dict['inner membrane'] = False
    if g[-1] in ['OM', 'LPO']:
        _dict['outer membrane'] = True
    else:
        _dict['outer membrane'] = False

    if (len(data[data['gene_name'] == g[0]]) == 0) & (len(data[data['b_number'] == g[1]]) == 0):
        continue

    for k, v in _dict.items():
        data.loc[data['gene_name'] == g[0], k] = v
        data.loc[data['b_number'] == g[1], k] = v

data.loc[data['cog_letter'].isin(
    ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q']), 'metabolism'] = True
data.loc[data['gene_name'].str.lower().isin(ribo_prots), 'ribosomal'] = True
data.to_csv('../../../data/literature/collated_mass_fractions.csv', index=False)

# %%
# Create aggregate summaries
envelope = data[data['envelope'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
envelope['localization'] = 'envelope'

membrane = data[data['membrane'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
membrane['localization'] = 'membrane'

periplasm = data[data['periplasm'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
periplasm['localization'] = 'periplasm'

inner_mem = data[data['inner membrane'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
inner_mem['localization'] = 'inner membrane'

outer_mem = data[data['outer membrane'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
outer_mem['localization'] = 'outer membrane'

cytoplasm = data[data['envelope'] == False].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
cytoplasm['localization'] = 'cytoplasm'

ribosomal = data[data['ribosomal'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
ribosomal['localization'] = 'ribosomal sector'

metabolic = data[data['metabolism'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
metabolic['localization'] = 'metabolic sector'

aggregated = pd.concat([envelope, membrane, periplasm, inner_mem, outer_mem,
                        cytoplasm, ribosomal, metabolic, inner_mem, outer_mem], sort=False)
aggregated.to_csv('../../../data/literature/summarized_mass_fractions.csv')
