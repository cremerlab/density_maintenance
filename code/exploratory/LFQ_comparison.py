#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from tqdm import tqdm
import size.viz 
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('./LFQ_wt.txt', delimiter = '\t')
data.head()
data['gene_names'] = [k.split('|')[-1] for k in data['Majority protein IDs']]
data = data[['Intensity 80', 'Intensity 99', 'Intensity 82', 'Intensity 83', 
             'LFQ intensity 80', 'LFQ intensity 99', 'LFQ intensity 82', 
             'LFQ intensity 83', 'gene_names']]

# Load the uniprot identifiers
uniprot = pd.read_csv('./uniprot_identifiers.csv')
mapper = {g[0]:g[1:] for g, _ in uniprot.groupby(['Entry name', 
                                                  'Gene names  (primary )',
                                                  'Gene names'])}
# Map
for k, v in tqdm(mapper.items()):
    data.loc[data['gene_names'] == k, 'primary_name'] = v[0]
    data.loc[data['gene_names'] == k, 'synonyms'] = v[1]

# Drop unmapped genes
data.dropna(inplace=True)

cond_mapper = {'80': ['glucose', 1], '99': ['glucose', 2], 
               '82': ['glucose+acetate', 1], '83': ['glucose+acetate', 2]}
dfs = pd.DataFrame([])
for k, v in cond_mapper.items():
    _df = data[['primary_name', 'synonyms',f'Intensity {k}', f'LFQ intensity {k}']].copy()
    _df['carbon_source'] = v[0]
    _df['replicate'] = v[1]
    _df.rename(columns={f'Intensity {k}': 'intensity',
                       f'LFQ intensity {k}': 'lfq_intensity'}, inplace=True)
    _df['intensity_frac'] = _df['intensity']/_df['intensity'].sum()
    _df['lfq_intensity_frac'] = _df['lfq_intensity']/_df['lfq_intensity'].sum()
    dfs = pd.concat([dfs, _df]) 

# Add categorization
loc_mapper = pd.read_csv('../../data/literature/Babu2018/Babu2018_full_classification.csv')
rename = {'IM': 'inner membrane',
          'LPI': 'inner membrane',
          'OM': 'outer membrane',
          'LPO': 'outer membrane',  
          'EC': 'extracellular',
          'MR': 'membrane related',
          'PE': 'periplasm',
          'CP': 'cytoplasm'}
for g, _ in tqdm(loc_mapper.groupby(['gene', 'location'])):
    dfs.loc[dfs['primary_name'].str.lower()==g[0].lower(), 'localization'] = rename[g[1]]

# Drop unmapped        
dfs.dropna(inplace=True)

ribo_prots =  ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg',
'rpsh', 'rpsi', 'rpsj', 'rpsk', 'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq',
'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple',
'rplf', 'rplj', 'rpll', 'rpli', 'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq',
'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma', 'rpmb',
'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

dfs.loc[dfs['primary_name'].str.lower().isin(ribo_prots), 'localization'] = 'ribosome'

#%%
agged = dfs.groupby(['carbon_source', 
                     'replicate',
                     'localization'])[['intensity_frac',
                               'lfq_intensity_frac']].sum(numeric_only=True
                                                          ).reset_index()

#%%
ms_data = pd.read_csv('../../data/collated/collated_mass_fractions.csv')
ms_data.loc[ms_data['gene_name'].str.lower().isin(ribo_prots), 'localization'] = 'ribosome'
#%%
fig, ax = plt.subplots(1, 3, figsize = (6, 2))

labeled = False
for g, d in ms_data.groupby(['source', 'condition', 'growth_rate_hr']):
    rib = d[d['localization'] == 'ribosome']
    mem = d[d['localization'].isin(['inner membrane', 'outer membrane', 'membrane related'])]
    peri = d[d['localization'] == 'periplasm']

    fmt = size.viz.style_point(g[0])    
    if g[0] == 'This Study':
        fmt['label'] = 'TMT Intensity'
    else:
        fmt['label'] = '__nolegend__'
    if labeled == True:
        fmt['label'] = '__nolegend__'
    
    ax[0].plot(rib['growth_rate_hr'].values[0], 
               rib['mass_frac'].sum(), **fmt)
    ax[1].plot(mem['growth_rate_hr'].values[0], 
               mem['mass_frac'].sum(), **fmt)
    ax[2].plot(peri['growth_rate_hr'].values[0], 
               peri['mass_frac'].sum(), **fmt)
    if g[0] == 'This Study':
        labeled = True 
# Map the growth rates
lam = ms_data[(ms_data['source'] == 'This Study') & 
(ms_data['condition'].isin(['glucose',
'glucose+acetate']))].groupby(['condition', 'growth_rate_hr',
'replicate'])['mass_frac'].sum().reset_index()

labeled = False
for g, d in agged.groupby(['carbon_source', 'replicate']):
    _lam = lam[(lam['condition'] == g[0]) & (lam['replicate'] == g[1])]['growth_rate_hr'].values[0]
    rib = d[d['localization'] == 'ribosome']
    mem = d[d['localization'].isin(['inner membrane', 'outer membrane', 'membrane related'])]
    peri = d[d['localization'] == 'periplasm']
    fmt = size.viz.style_point('This Study')
    fmt['marker'] = '*'
    fmt['color'] = cor['primary_red']
    fmt['markersize'] = 6
    if labeled == False:
        fmt['label'] = 'LFQ-intensity' 
    else:
        fmt['label'] = '__nolegend__'
    ax[0].plot(_lam, rib['lfq_intensity_frac'].sum(), **fmt)
    ax[1].plot(_lam, mem['lfq_intensity_frac'].sum(), **fmt)
    ax[2].plot(_lam, peri['lfq_intensity_frac'].sum(), **fmt)
    fmt['marker'] = 'v'
    fmt['color'] = cor['primary_blue']
    fmt['markersize'] = 4
    if labeled == False:
        fmt['label'] = 'Intensity'
    else:
        fmt['label'] = '__nolegend__'
    ax[0].plot(_lam, rib['intensity_frac'].sum(), **fmt)
    ax[1].plot(_lam, mem['intensity_frac'].sum(), **fmt)
    ax[2].plot(_lam, peri['intensity_frac'].sum(), **fmt)
    labeled = True


for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].legend(fontsize=6, loc='upper left')
ax[0].set_ylabel('ribosomal mass fraction', fontsize=6)
ax[1].set_ylabel('membrane mass fraction', fontsize=6)
ax[2].set_ylabel('periplasmic mass fraction', fontsize=6)












