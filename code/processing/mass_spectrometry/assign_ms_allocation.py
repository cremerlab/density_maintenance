#%%
import pandas as pd 

# Load the mass spectrometry data from first data collection run.
ms_data = pd.read_csv('processed_mass_fractions.csv') 

# Drop invalid lacZ measurements
ms_data = ms_data[ms_data['strain']!='lacZ']

# Load and collate the mass spectrometry data from the second collection run. 
ms_data = pd.concat([ms_data, pd.read_csv('./processed_mass_fractions_ppGpp_glucoseCAA.csv')])

#%%
# Load the gene classification and define the ribosomal proteins. 
gene_class = pd.read_csv('../../../data/literature/genes_classification_all.csv')

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
    d['cog_category'] = lcz['COG_functionname'].values[0]
    d['cog_letter'] = lcz['COG_function'].values[0]
    filt = pd.concat([filt, d])

filt = filt[['strain', 'carbon_source', 'date', 'inducer_conc', 'replicate',
'name', 'mass_frac', 'localization', 'cog_category', 'cog_letter']]
filt.to_csv('./mass_spectrometry_localization.csv', index=False) 

#%%
# Do two passes of computing the allocation. First based on localization, second 
# on ribosomal content
groups = ['strain', 'carbon_source', 'date', 'inducer_conc',
          'replicate', 'localization']
allocs = filt.groupby(groups)['mass_frac'].sum().reset_index()
phi_rib = filt[filt['name'].str.lower().isin(ribo_prots)].groupby(groups[:-1])['mass_frac'].sum().reset_index()
phi_rib['localization'] = 'phi_rib'

# Concatenate and then pivot
concat = pd.concat([allocs, phi_rib])
concat.rename(columns={'localization':'allocation'}, inplace=True)
concat.to_csv('./mass_spectrometry_allocation_long.csv', index=False)

#%%
groups[-1] = 'allocation'
groups.append('mass_frac')
pivot = pd.DataFrame([])
for g, d in concat.groupby(groups[:-2]):     
    _df = pd.DataFrame([d.values[0]], columns=groups)
    for _g, _d in d.groupby('allocation'):
        _df[_g] = _d['mass_frac'].values 
    _df.drop(columns=['allocation', 'mass_frac'], inplace=True)
    pivot = pd.concat([pivot, _df])
pivot.to_csv('./mass_spectrometry_allocation_wide.csv', index=False)
