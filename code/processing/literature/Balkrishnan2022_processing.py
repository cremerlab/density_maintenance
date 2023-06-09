# %%
import numpy as np
import pandas as pd
raw = pd.read_csv(
    '../../../data/literature/Balakrishnan2022/Balakrishnan2022_rna_fraction_raw.csv')
mapper = pd.read_csv(
    '../../../data/literature/Balakrishnan2022/Balakrishnan2022_condition_mapper.csv')


# Standardize gene names
gene_list = pd.read_csv(
    '../../../data/literature/Belliveau2021/ecoli_genelist_master.csv')
gene_dict = {}
for g, d in gene_list.groupby(['b_number',
                               'gene_name',
                               'go_terms',
                               'cog_letter',
                               'cog_class',
                               'cog_desc']):
    gene_dict[g[1].lower()] = {'b_number': g[0],
                               'go_terms': g[2],
                               'cog_letter': g[3],
                               'cog_class': g[4],
                               'cog_desc': g[5]}
    gene_dict[g[0].lower()] = {'gene_name': g[1],
                               'go_terms': g[2],
                               'cog_letter': g[3],
                               'cog_class': g[4],
                               'cog_desc': g[5]}
names = []
go_terms = []
cog_letter = []
cog_class = []
unclassed = 0
unclassed_genes = []
unclassed_bnumbers = []
for g, b in zip(raw['gene'].values, raw['locus'].values):
    try:
        _info = gene_dict[b]
        if len(_info['gene_name']) > 5:
            names.append(g)
        else:
            _g = _info['gene_name']
            if len(_g) > 3:
                _g[3:] = _g[3:].upper()
            names.append(_g['gene_name'])
        cog_letter.append(_info['cog_letter'])
        cog_class.append(_info['cog_class'])
        go_terms.append(_info['go_terms'])
    except:
        try:
            _info = gene_dict[g.split('_')[0].lower()]
            names.append(g.split('_')[0])
            cog_letter.append(_info['cog_letter'])
            cog_class.append(_info['cog_class'])
            go_terms.append(_info['go_terms'])
        except:
            print(f'Could not map {g}:{b}')
            unclassed += 1
            unclassed_genes.append(g)
            unclassed_bnumbers.append(b)
            for ell in [names, go_terms, cog_letter, cog_class]:
                ell.append(np.nan)

raw['gene_name'] = names
raw['go_terms'] = go_terms
raw['cog_letter'] = cog_letter
raw['cog_class'] = cog_class
raw.dropna(inplace=True)
raw.drop(columns=['gene length (nt)', 'gene'], inplace=True)
melted = raw.melt(id_vars=['gene_name', 'go_terms',
                  'cog_class', 'cog_letter', 'locus',])
melted.rename(columns={'variable': 'condition',
              'locus': 'b_number',
                       'value': 'mass_frac'}, inplace=True)
melted['dataset_name'] = 'Balakrishnan et al. 2022'

for g, _ in mapper.groupby(['Sample ID', 'Strain',  'Growth rate (1/h)', 'Carbon source', 'replicate', 'Supplement']):
    if g[0][0] == 'c':
        csource = g[3] + f' (c limitation + {g[-1]})'  # , replicate {g[-2]}) '
    elif g[0][0] == 'r':
        csource = g[3] + f' (r limitation + {g[-1]})'  # , replicate {g[-2]})'
    elif g[0][0] == 'a':
        csource = g[3] + f' (n limitation + {g[-1]})'  # , replicate {g[-2]})'
    melted.loc[melted['condition'] == g[0], 'strain'] = g[1]
    melted.loc[melted['condition'] == g[0], 'growth_rate_hr'] = g[2]
    melted.loc[melted['condition'] == g[0], 'condition'] = csource
# melted.drop(columns=['condition'], inplace=True)
melted = melted[melted['condition'].str.contains('(c limitation)')]
print(len(melted))
melted = melted.groupby(['gene_name', 'go_terms', 'cog_class', 'cog_letter',
                        'b_number', 'condition', 'dataset_name', 'strain']).mean().reset_index()
print(len(melted))
melted.to_csv(
    '../../../data/literature/Balakrishnan2022/Balakrishnan2022_processed.csv', index=False)
