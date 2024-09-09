#%%
import numpy as np 
import pandas as pd 

# Define the valid samples to use (namely, C-limitation)
table8_valid_idx = ['Lib-06',
                    'Lib-24',
                    'Lib-25',
                    'Lib-26',
                    'Lib-27',
                    'Lib-28',
                    'Lib-29',
                    'Lib-30']

# Define the valid samples to use (namely, C-limitation)
table9_valid_idx = ['A1-1',
                    'A1-2',
                    'A1-3',
                    'F1-1'
                    'C2'
                    'C3',
                    'C4',
                    'C5',
                    'C6',
                    'C7',
                    'C8',
                    'D6',
                    'D7',
                    'D8',
                    'F4',
                    'F5',
                    'F6',
                    'F7', 
                    'F8']

# Load the sample idx csv to make a LUT
table8_keys = pd.read_csv('../../../data/literature/Mori2021/Mori2021_table8_idx.csv')
table9_keys = pd.read_csv('../../../data/literature/Mori2021/Mori2021_table9_idx.csv')

#%%
table8_dict = {} 
for g, d in table8_keys.groupby(['Sample ID', 'Carbon Source', 
                                 'Strain', 'Doubling time (min)']):
    table8_dict[g[0]] = {'condition': g[1],
                      'strain': g[2],
                      'growth_rate_hr': np.log(2) / (g[-1] / 60)}

table9_dict = {}
for g, d in table9_keys.groupby(['Sample ID', 'Strain', 'Growth rate (1/h)', 'Carbon source']):
    table9_dict[g[0]] = {'condition': g[3],
                         'strain': g[1],
                         'growth_rate_hr':g[2]}

# Load and melt the tables
table8 = pd.read_csv('../../../data/literature/Mori2021/Mori2021_table8.csv')
table9 = pd.read_csv('../../../data/literature/Mori2021/Mori2021_table9.csv')
table8_melted = table8.melt(['Gene name', 'Gene locus', 'Protein ID'], var_name='idx', value_name='mass_frac')
table9_melted = table9.melt(['Gene name', 'Gene locus', 'Protein ID'], var_name='idx', value_name='mass_frac')

# Restrict the datatable to only those valid idx
table8_valid = table8_melted[table8_melted['idx'].isin(table8_valid_idx)]
table9_valid = table9_melted[table9_melted['idx'].isin(table9_valid_idx)]

# Populate the tables with identifiers captured in the dict
for tab, keys in zip([table8_valid, table9_valid], [table8_dict, table9_dict]):
    for val in ['condition', 'strain', 'growth_rate_hr']:
        _vals = [keys[idx][val] for idx in tab['idx'].values]
        tab[val] = _vals

# Merge the fixed tables
merged = pd.concat([table8_valid, table9_valid])

# Rename the columns to the standard 
merged = merged[['Gene name', 'Gene locus', 'mass_frac', 'condition', 'strain', 'growth_rate_hr']]
merged.rename(columns={'Gene name': 'gene_name', 'Gene locus':'b_number'}, inplace=True)
merged['dataset_name'] = 'Mori et al. 2021'
merged = merged[merged['mass_frac'] > 0]

# %%
# Load the master gene list and make a LUT
gene_list = pd.read_csv('../../../data/literature/Belliveau2021/ecoli_genelist_master.csv')
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

    gene_dict[g[0]] = {'gene_name': g[1],
                       'go_terms':g[2],
                       'cog_letter': g[3],
                       'cog_class': g[4],
                       'cog_desc': g[5]}


# Load the master gene list for consistent annotation
unmapped = []
for v in ['go_terms', 'cog_class', 'cog_letter', 'cog_desc']:
    info = []
    for i, b in enumerate(merged['gene_name'].values):
        try:
            info.append(gene_dict[b.lower().split('_')[0]][v])
        except KeyError:
            info.append(gene_dict[merged['b_number'].values[i]][v])
    merged[v] = info

# Deal with the case where there are replicates by computing the mean for each gene
merged = merged.groupby(['gene_name', 'b_number', 'condition', 'strain', 'growth_rate_hr',
'dataset_name', 'go_terms', 'cog_class', 'cog_letter', 'cog_desc']).mean().reset_index()
merged
# %%
merged.to_csv('../../../data/literature/Mori2021/Mori2021_processed.csv', index=False)


