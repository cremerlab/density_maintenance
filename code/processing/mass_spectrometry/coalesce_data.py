#%%
import numpy as np 
import pandas as pd 

gene_class = pd.read_csv('../../../data/literature/genes_classification_all.csv')
size_data = pd.read_csv('./aggregated_size_measurements.csv')
ms_data = pd.read_csv('./processed_mass_fractions.csv')
ms_data = pd.concat([ms_data, pd.read_csv('./processed_mass_fractions_ppGpp_glucoseCAA.csv')])
growth_data = pd.read_csv('./collated_growth_rates.csv')
# RP_data = pd.read_csv('./collated_protein_RNA_measurements.csv')
# RP_data.rename(columns={'date_collected': 'date'}, inplace=True)
# RP_data['RNA_to_protein'] = RP_data['ug_RNA_per_od_600nm'] / RP_data['ug_protein_per_od_600nm']


#%%
# Define the ribosomal proteins
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

# Define the localization
locs = {'membrane': ['IM', 'LPI', 'LPO', 'OM', 'MR'],
       'periplasm': ['PE', 'EC'],
       'cytoplasm': ['CP'],       
}
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
    elif lcz['location'].values[0] in locs['cytoplasm']:
        d['localization'] = 'phi_cyto'
    else:
        d['localization'] = 'unassigned_localization'
    if g[0].lower() in ribo_prots:
        d['sector'] = 'phi_rib'
    elif g[0] in ['lacZ']:
        d['sector'] = 'phiX-lacZ'
    elif g[0] in ['relA']:
        d['sector'] = 'phiX-relA'
    else:
        d['sector'] = 'unassigned_sector'
    d['cog_category'] = lcz['COG_functionname'].values[0]
    d['cog_letter'] = lcz['COG_function'].values[0]
    filt = pd.concat([filt, d])
filt.to_csv('./mass_spec_sector_assignments.csv', index=False) 

#%%
agged = pd.DataFrame([])
groupby = ['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate']
for g, d in filt.groupby(groupby):
    grp = d.groupby(['localization', 'sector'])['mass_frac'].sum().reset_index()
    for i, v in enumerate(['localization', 'sector']):
        _grp = grp.groupby(v)['mass_frac'].sum().reset_index()
        if i == 0:
            _df = pd.DataFrame([{j:k for j, k in zip(_grp[v].values, _grp['mass_frac'].values)}])
        else:
            for j, k in zip(_grp[v].values, _grp['mass_frac'].values):
                _df[j] = k
    for g, k in zip(groupby, g): 
        _df[g] = k
    agged = pd.concat([agged, _df])
#%%
#%%
keys = ['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate', 'localization', 'sector']
# agged_ms = filt.groupby(keys)[['mass_frac', 'rel_mass_frac']].sum().reset_index()

# agged_ms_rel = agged_ms.pivot(index=keys[:-1], columns=['localization', 'sector'], values='rel_mass_frac').reset_index()
# agged_ms_rel.to_csv('./mass_spec_aggregated_sectors_relative.csv', index=False)

# agged_ms = agged_ms.pivot(index=keys[:-1], columns=['localization', 'sector'], values='mass_frac').reset_index()
agged.to_csv('./mass_spec_aggregated_sectors.csv', index=False)

#%%
merged = agged.merge(size_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged =  merged.merge(growth_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged = merged.merge(RP_data, on=['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged.to_csv(f'./total_collated_data.csv', index=False)


