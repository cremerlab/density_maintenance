#%%
import numpy as np

import pandas as pd

# Load in the total RNA + protein data 
rna_prot = pd.read_csv('../RNA_protein_quantification/rna_protein_per_biomass.csv')

# Load in the cells per biomass counts 
counts = pd.read_csv('../flow_cytometry/cells_per_biomass.csv')

# Compute the mean and std for the cells per biomass
agged_counts = counts.groupby(['carbon_source', 'growth_rate_hr'])['cells_per_biomass'].agg(('mean', 'std')).reset_index()
mapper = {g[0]:g[1] for g, _ in agged_counts.groupby(['carbon_source', 'mean'])}
for k, v in mapper.items():
    rna_prot.loc[rna_prot['carbon_source'] == k, 'mean_cells_per_biomass'] = v

rna_prot['fg_rna_per_cell'] = rna_prot['ug_rna_per_biomass'] * 1E9 / rna_prot['mean_cells_per_biomass']
rna_prot['fg_protein_per_cell'] = rna_prot['ug_protein_per_biomass'] * 1E9 / rna_prot['mean_cells_per_biomass']
rna_prot = rna_prot[['strain', 'carbon_source', 'growth_rate_hr', 'fg_rna_per_cell', 'fg_protein_per_cell']]
rna_prot.to_csv('../../../data/collated/experimental_rna_protein_per_cell.csv', index=False)
