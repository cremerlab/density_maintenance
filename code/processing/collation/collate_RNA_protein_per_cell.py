#%%
import numpy as np
import pandas as pd

# Load in the growth rate measurements from the mass spectrometry samples.
growth_rates = pd.read_csv('../mass_spectrometry/compiled_data/aggregated_growth_measurements.csv')
growth_rates = growth_rates[growth_rates['strain']=='wildtype']
growth_rates = growth_rates.groupby('carbon_source')['growth_rate_hr'].mean().reset_index()
 
# Create a hashmap
mapper = {g[0]:g[1] for g, _ in growth_rates.groupby(['carbon_source', 'growth_rate_hr'])}


# Load in the total RNA + protein data 
rna_prot = pd.read_csv('../RNA_protein_quantification/rna_protein_per_biomass.csv')

# Load in the cells per biomass counts 
counts = pd.read_csv('../flow_cytometry/cells_per_biomass.csv')

# Add in the growth rates to the datasets. 
for k, v in mapper.items():
    for d in [rna_prot, counts]:
        d.loc[d['carbon_source']==k, 'growth_rate_hr'] = v

# Save all to disk before merger and aggregation
rna_prot.to_csv('../../../data/collated/experimental_rna_protein_data.csv', index=False)
counts.to_csv('../../../data/collated/experimental_cells_per_biomass_data.csv', index=False)

# Compute the mean and std for the cells per biomass
agged_counts = counts.groupby(['carbon_source', 'growth_rate_hr'])['cells_per_biomass'].agg(('mean', 'std')).reset_index()

# Compute the mean and std for each quantity separately 
agged_rna = rna_prot.groupby(['carbon_source', 'growth_rate_hr'])['ug_rna_per_biomass'].agg(('mean', 'std')).reset_index()
agged_prot = rna_prot.groupby(['carbon_source', 'growth_rate_hr'])['ug_protein_per_biomass'].agg(('mean', 'std')).reset_index()

# Include the new means
agged_rna['fg_rna_per_cell_mu'] = 1E9 * agged_rna['mean'] / agged_counts['mean']
agged_prot['fg_protein_per_cell_mu'] = 1E9 * agged_prot['mean'] / agged_counts['mean']

# Propagate the errors
agged_rna['fg_rna_per_cell_sigma'] = np.sqrt((1E9 * agged_rna['std']/agged_counts['mean'])**2 +\
                                              (1E9 * agged_rna['mean'] * agged_counts['std'])**2/(agged_counts['mean']**4))
agged_rna.drop(columns=['mean', 'std'], inplace=True)
agged_prot['fg_protein_per_cell_sigma'] = np.sqrt((1E9 * agged_prot['std']/agged_counts['mean'])**2 +\
                                              (1E9 * agged_prot['mean'] * agged_counts['std'])**2/(agged_counts['mean']**4))
agged_prot.drop(columns=['mean', 'std'], inplace=True)

# Merge and save
merged  = agged_rna.merge(agged_prot, on=['carbon_source', 'growth_rate_hr'], how='inner')
merged.to_csv('../../../data/collated/experimental_rna_protein_per_cell.csv', index=False)


