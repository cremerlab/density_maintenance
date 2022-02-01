#%%
import numpy as np 
import pandas as pd 

# Load the data and dropped unnecessary columns
data = pd.read_csv('../../data/source/Soufi2015/soufi2015_longform_annotated.csv')
data.drop(columns=['Unnamed: 0', 'reported_tot_per_cell', 'corrected_volume', 'dataset'],
                inplace=True)

# Compute the mass fractions
dfs = []
for g, d in data.groupby(['growth_rate_hr']):
    tot = d['reported_fg_per_cell'].sum()
    d['mass_frac'] = d['reported_fg_per_cell'].values / tot
    dfs.append(d)

# Concatenate the dataset and drop an unnecessary column
data = pd.concat(dfs, sort=False)
data.drop(columns=['reported_fg_per_cell'], inplace=True)

# Load the master gene list and make a LUT
gene_list = pd.read_csv('../../data/source/Belliveau2021/ecoli_genelist_master.csv')
gene_dict = {}

for g, d in gene_list.groupby(['b_number', 
                               'gene_name', 
                               'go_terms', 
                               'cog_letter', 
                               'cog_class',
                               'cog_desc']):
    gene_dict[g[1].lower()] = {
                       'go_terms':g[2],
                       'cog_letter': g[3],
                       'cog_class': g[4],
                       'cog_desc': g[5]}

# Add the stuff in the gene dict to the Soufi dataframe
genes = data['gene_name'].values
for k in ['go_terms', 'cog_letter', 'cog_class', 'cog_desc']:
    data[k] = [gene_dict[g.lower().split('_')[0]][k] for g in genes]

# Save
data.to_csv('../../data/source/Soufi2015/Soufi2015_processed.csv', index=False)
# %%
