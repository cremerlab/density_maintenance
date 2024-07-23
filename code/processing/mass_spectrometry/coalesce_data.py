#%%
import numpy as np 
import pandas as pd 

gene_class = pd.read_csv('../../../data/literature/genes_classification_all.csv')
size_data = pd.read_csv('./aggregated_size_measurements.csv')
ms_data = pd.read_csv('./processed_mass_fractions.csv')
growth_data = pd.read_csv('./collated_growth_rates.csv')
RP_data = pd.read_csv('./collated_protein_RNA_measurements.csv')
RP_data.rename(columns={'date_collected': 'date'}, inplace=True)

# Define the ribosomal proteins
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

# Define the localization
locs = {'membrane': ['IM', 'LPI', 'LPO', 'OM'],
       'periplasm': ['EC', 'PE']}
filt = pd.DataFrame([])
for g, d in ms_data.groupby(['name', 'synonyms']):
    lcz = gene_class.loc[(gene_class['gene']==g[0]) | (gene_class['gene'].isin(g[1].split()))]
    if len(lcz) == 0:
        print(f'could not classify {g}')     
        continue
    if lcz['location'].values[0] in locs['membrane']:
        d['localization'] = 'phi_mem'
    elif lcz['location'].values[0] in locs['periplasm']:
        d['localization'] = 'phi_peri'
    elif g[0].lower() in ribo_prots:
        d['localization'] = 'phi_rib'
    elif g[0] in ['lacZ']:
        d['localization'] = 'phiX-lacZ'
    elif g[0] in ['relA']:
        d['localization'] = 'phiX-relA'
    else:
        continue 
    filt = pd.concat([filt, d])
filt.to_csv('./mass_spec_sector_assignments.csv', index=False) 

#%%
keys = ['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate', 'localization']
agged_ms = filt.groupby(keys)['mass_frac'].sum().reset_index()
agged_ms = agged_ms.pivot(index=keys[:-1], columns='localization', values='mass_frac').reset_index()
agged_ms.to_csv('./mass_spec_aggregated_sectors.csv', index=False)

#%%
merged = agged_ms.merge(size_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged =  merged.merge(growth_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged = merged.merge(RP_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged.to_csv('./total_collated_data.csv', index=False)


