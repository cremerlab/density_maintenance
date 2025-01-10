#%%
import pandas as pd 

# Load and colalted literature data
files = ['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']  
dfs = [pd.read_csv(
    f'../../../data/literature/{f}/{f}_processed.csv') for f in files]
lit_data = pd.concat(dfs, sort=False)

# Restrict and simplify the columns
lit_data = lit_data[['gene_name', 'dataset_name', 'strain', 'condition', 
             'growth_rate_hr', 'cog_class', 'cog_letter', 'mass_frac']]
lit_data['replicate'] = 0
lit_data.rename(columns={'dataset_name': 'source'}, inplace=True)
lit_data['gene_name'] = [g.lower() for g in lit_data['gene_name'].values]

#%%
# Load the gene classification and define the ribosomal proteins. 
gene_class = pd.read_csv('../../../data/literature/genes_classification_all.csv')
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg',
'rpsh', 'rpsi', 'rpsj', 'rpsk', 'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq',
'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple',
'rplf', 'rplj', 'rpll', 'rpli', 'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq',
'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma', 'rpmb',
'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

# Define the localization
locs = {'membrane': ['IM', 'LPI', 'LPO', 'OM', 'MR'],
       'periplasm': ['PE', 'EC'],
       'cytoplasm': ['CP'],       
}
filt = pd.DataFrame([])
# Count which are not classified
not_classified = []
for g, d in lit_data.groupby('gene_name'):
    lcz = gene_class.loc[gene_class['gene'].str.lower()==g.lower()]
    if len(lcz) == 0:
        print(f'could not classify {g}')     
        not_classified.append(g)
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

filt = filt[['gene_name', 'source', 'condition', 'growth_rate_hr', 'mass_frac', 'localization', 'cog_letter', 'cog_class', 'strain']]
filt.rename(columns={'gene_name':'name'}, inplace=True)
filt.to_csv('./compiled_data/compiled_literature_mass_fractions.csv', index=False)

#%%
# Do two passes of computing the allocation. First based on localization, second 
# on ribosomal content
groups = ['source', 'strain', 'condition', 'growth_rate_hr', 'localization']
allocs = filt.groupby(groups)['mass_frac'].sum().reset_index()
phi_rib = filt[filt['name'].str.lower().isin(ribo_prots)].groupby(groups)['mass_frac'].sum().reset_index()
phi_rib['localization'] = 'phi_rib'

# Concatenate and then pivot
concat = pd.concat([allocs, phi_rib])
concat.rename(columns={'localization':'allocation'}, inplace=True)
concat.to_csv('./compiled_data/compiled_literature_allocation_assigments_long.csv', index=False)

# Generate a wide-format table
groups[-1] = 'allocation'
groups.append('mass_frac')
pivot = pd.DataFrame([])
for g, d in concat.groupby(groups[:-2]):     
    _df = pd.DataFrame([d.values[0]], columns=groups)
    for _g, _d in d.groupby('allocation'):
        _df[_g] = _d['mass_frac'].values 
    _df.drop(columns=['allocation', 'mass_frac'], inplace=True)
    pivot = pd.concat([pivot, _df])
pivot.to_csv('../../../data/collated/compiled_literature_allocation_assigments_wide.csv', index=False)