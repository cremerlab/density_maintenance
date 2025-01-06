#%%
import pandas as pd 

# Load the mass spec, size, and growth measurements together
size_data = pd.read_csv('../aggregated_size_measurements.csv')
growth_data = pd.read_csv('../aggregated_growth_measurements.csv')
alloc = pd.read_csv('../mass_spectrometry_allocation_wide.csv')

# Merge on the shared columns
shared = ['strain', 'carbon_source', 'inducer_conc', 'replicate', 'date']
merged = size_data.merge(growth_data, on=shared, how='inner')
merged = merged.merge(alloc, on=shared, how='inner')
merged.to_csv('../total_aggregated_data.csv')


# Set a mapper from strain, date, inducer_conc, replicate to growth rate. 
groups = ['strain', 'inducer_conc', 'date', 'replicate', 'growth_rate_hr']
mapper = {g[:-1]:g[-1] for g, _ in growth_data.groupby(groups)}

# Add growth rate information to the raw mass spec data. 
ms_data = pd.read_csv('../mass_spectrometry_localization.csv')
for k, v in mapper.items():
    ms_data.loc[(ms_data[groups[0]]==k[0]) & 
                (ms_data[groups[1]]==k[1]) & 
                (ms_data[groups[2]]==k[2]) & 
                (ms_data[groups[3]]==k[3]), 'growth_rate_hr'] = v
ms_data.to_csv('../mass_spectrometry_mass_fractions.csv', index=False)

