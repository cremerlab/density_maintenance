#%%
import numpy as np 
import pandas as pd 
import tqdm

# Load and colalted literature data
files = ['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']  
dfs = [pd.read_csv(
    f'../../../data/literature/{f}/{f}_processed.csv') for f in files]
lit_data = pd.concat(dfs, sort=False)

# Exclude valgepea et al. which is not well designed to measure membrane proteins
lit_data = lit_data[lit_data['dataset_name'] != 'Valgepea et al. 2013']

# Restrict and simplify the columns
lit_data = lit_data[['gene_name', 'dataset_name', 'strain', 'condition', 
             'growth_rate_hr', 'cog_class', 'cog_letter', 'mass_frac']]
lit_data['replicate'] = 0
lit_data.rename(columns={'dataset_name': 'source'}, inplace=True)

#%% 
# From the literature data, make map between gene names and cog info
cog_mapper = {g[0].lower():[g[1], g[2]] for g, _ in lit_data.groupby(['gene_name', 'cog_class', 'cog_letter'])}

#%%
# Load our data and restrict
ms_data = pd.read_csv('../mass_spectrometry/processed_mass_fractions.csv')
lam_data = pd.read_csv('../mass_spectrometry/collated_growth_rates.csv')
 
merged_data = pd.merge(ms_data, lam_data, on=['date', 'strain', 'carbon_source',
                                              'inducer_conc', 'replicate'], how='outer')
merged_data = merged_data[merged_data['strain']=='wildtype']

# Make columns consistent
merged_data['strain'] = 'NCM3722'
merged_data['condition'] = merged_data['carbon_source']
merged_data['source'] = 'This Study'
merged_data.rename(columns={'name':'gene_name'}, inplace=True)
merged_data = merged_data[['gene_name', 'source', 'strain', 'condition', 'growth_rate_hr',
                           'mass_frac', 'replicate']] 
# Add cog info to our data
unmapped = 0
for g, d in merged_data.groupby('gene_name'):
    if g.lower() not in cog_mapper.keys():
        print(f'{g} not in cog_mapper')
        unmapped += 1
        continue
    cog_class, cog_letter = cog_mapper[g.lower()]
    merged_data.loc[merged_data['gene_name']==g, 'cog_class'] = cog_class
    merged_data.loc[merged_data['gene_name']==g, 'cog_letter'] = cog_letter
# Drop the unmapped genes
merged_data.dropna(inplace=True)

#%% 
# Link the data sets into a single dataframe
data = pd.concat([lit_data, merged_data], sort=False)

#%% Add localization
mapper = pd.read_csv('../../../data/literature/Babu2018/Babu2018_full_classification.csv')
rename = {'IM': 'inner membrane',
          'LPI': 'inner membrane',
          'OM': 'outer membrane',
          'LPO': 'outer membrane',  
          'EC': 'extracellular',
          'MR': 'membrane related',
          'PE': 'periplasm',
          'CP': 'cytoplasm'}
for g, _ in tqdm.tqdm(mapper.groupby(['gene', 'location'])):
    data.loc[data['gene_name'].str.lower()==g[0].lower(), 'localization'] = rename[g[1]]

# Save the data
# Keep only the 37° C data
data = data[data['condition'] != '42C']
data = data[~data['mass_frac'].isnull()]
data.to_csv('../../../data/collated/collated_mass_fractions.csv', index=False)


#%%
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']
# Aggregate and the literature data by type and localization.
agged = pd.DataFrame([])
for g, d in data.groupby(['source', 'condition', 'growth_rate_hr', 'strain', 'replicate']):
    d = d[d['mass_frac'] > 0]
    phi_cyt = d[d['localization']=='cytoplasm']['mass_frac'].values.sum()
    phi_rib = d[d['gene_name'].str.lower().isin(ribo_prots)]['mass_frac'].values.sum()
    phi_mem = d[d['localization'].isin(['inner membrane', 'outer membrane', 'membrane related'])]['mass_frac'].values.sum()
    phi_peri = d[d['localization']=='periplasm']['mass_frac'].values.sum()
    _df = pd.DataFrame({'source': g[0],
                        'condition':g[1],
                        'growth_rate_hr':g[2],
                        'strain':g[3],
                        'replicate': g[4],
                        'phi_cyt': phi_cyt,
                        'phi_rib': phi_rib,
                        'phi_mem': phi_mem,
                        'phi_peri': phi_peri},
                        index=[0])
    agged = pd.concat([agged, _df], sort=False)
agged.to_csv('../../../data/collated/literature_mass_spec_aggregated.csv', index=False)

