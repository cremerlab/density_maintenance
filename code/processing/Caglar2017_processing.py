#%%
import numpy as np 
import pandas as pd 

# LOad data and omit high salt condition
data = pd.read_csv('../../data/source/Caglar2017/caglar2017_longform_annotated.csv')
data = data[data['condition'] != 'NaCl_stress']

# Drop unnecessary columns
data = data[['gene_name', 'condition', 'growth_rate_hr', 'annotation', 
             'cog_category', 'cog_class', 'cog_letter', 'fg_per_cell', 
             'dataset', 'strain']]

# Compute the mass fraction
dfs = []
for g, d in data.groupby(['condition']):
    tot = d['fg_per_cell'].sum()
    d['mass_frac'] = d['fg_per_cell'].values / tot
    dfs.append(d)
data = pd.concat(dfs, sort=False)
# Drop fg per cell and rename dataset
data.drop(columns=['fg_per_cell'], inplace=True)
data.rename(columns={'dataset':'dataset_name'})

# %%
# Load the gene dict and generate a map to b number
gene_list = pd.read_csv('../../data/source/Belliveau2021/ecoli_genelist_master.csv')
gene_dict = {}

for g, d in gene_list.groupby(['b_number', 
                               'gene_name', 
                               'go_terms', 
                               'cog_letter', 
                               'cog_class',
                               'cog_desc']):
    gene_dict[g[1].lower()] = {'b_number': g[0],
                       'go_terms':g[2],
                       'cog_letter': g[3],
                       'cog_class': g[4],
                       'cog_desc': g[5]}
 
# %%
# Link the genes as necessary
data = data[data['mass_frac'] > 0]
for v in ['go_terms', 'cog_class', 'cog_letter', 'cog_desc']:
    info = []
    for i, b in enumerate(data['gene_name'].values):
        info.append(gene_dict[b.lower().split('_')[0]][v])
    data[v] = info
data['dataset_name'] = 'Caglar et al. 2017'
# %%
# Save
data.to_csv('../../data/source/Caglar2017/Caglar2017_processed.csv', index=False)


# %%
