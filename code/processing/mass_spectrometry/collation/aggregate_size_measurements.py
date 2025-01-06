#%%
import pandas as pd 
import glob

# Load all of the size measurements from the directory
size_data = pd.concat([pd.read_csv(f) for f in glob.glob('../2024-*/processed/*sizes.csv')])

# Compute the mean cell size quantities
groups = ['date', 'replicate', 'strain', 'carbon_source', 'inducer_conc']
quants = ['width_median', 'length', 'volume', 'surface_area', 'surface_to_volume']
agged = size_data.groupby(groups)[quants].mean().reset_index()

# Rename the columns appropriately
agged.rename(columns = {'width_median':'width_um',
                                  'length':'length_um',
                                  'volume': 'volume_fL',
                                  'surface_area': 'surface_area_um2',
                                  'surface_to_volume': 'surface_to_volume_inv_um'
                                  }, inplace=True)



# Load the file that prescribes valid experiments
metadata = pd.read_csv('../valid_experimental_metadata.csv')
metadata.rename(columns={'date_collected': 'date'}, inplace=True)

# Merge in only those experiments from the aggregated size measurements that are 
# deemed valid
groups = ['date', 'replicate', 'strain', 'carbon_source', 'inducer_conc']
valid = agged.merge(metadata, on=groups, how='inner')
valid.to_csv('../aggregated_size_measurements.csv', index=False)