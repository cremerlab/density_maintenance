#%%
import pandas as pd 
import glob 

# Load all of the size measurements from the directory
growth_data = pd.concat([pd.read_csv(f) for f in glob.glob('../2024-*/processed/*growth_rates.csv')])

#%%
# Compute the mean cell size quantities
groups = ['date', 'replicate', 'strain', 'carbon_source', 'inducer_conc']
quants = ['growth_rate_hr']
growth_data.loc[growth_data['strain'].isnull(), 'strain'] = 'wildtype'
growth_data.loc[growth_data['inducer_conc'].isnull(), 'inducer_conc'] = 0
agged = growth_data.groupby(groups)[quants].mean().reset_index()


# Load the file that prescribes valid experiments
metadata = pd.read_csv('../valid_experimental_metadata.csv')
metadata.rename(columns={'date_collected': 'date'}, inplace=True)

# Merge in only those experiments from the aggregated size measurements that are 
# deemed valid
groups = ['date', 'replicate', 'strain', 'carbon_source', 'inducer_conc']
valid = agged.merge(metadata, on=groups, how='inner')
valid.to_csv('../compiled_data/aggregated_growth_measurements.csv', index=False)