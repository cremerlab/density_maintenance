#%% 
import pandas as pd 

files =['Mori2021', 'Soufi2015', 'Caglar2017', 'Belliveau2021']
dfs = [pd.read_csv(f'../../data/source/{f}/{f}_processed.csv') for f in files]
data = pd.concat(dfs, sort=False)
data.drop(columns=['cog_category', 'cog_desc', 'gene_product', 'annotation', 'dataset'],
            inplace=True)
data.to_csv('../../data/compiled_mass_fractions.csv', index=False)
