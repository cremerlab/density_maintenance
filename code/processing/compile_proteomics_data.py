# %%
import pandas as pd

files = ['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']
dfs = [pd.read_csv(
    f'../../data/literature/{f}/{f}_processed.csv') for f in files]
data = pd.concat(dfs, sort=False)

# Categorize
data['periplasm'] = False
data.loc[(data['go_terms'].str.contains('GO:0042597'))
         | (data['go_terms'].str.contains('GO:0005620')), 'periplasm'] = True

# %%
data.drop(columns=['b_number', 'cog_category', 'cog_desc', 'gene_product', 'annotation', 'dataset'],
          inplace=True)
data.to_csv('../../data/compiled_mass_fractions.csv', index=False)

# %%
