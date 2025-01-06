#%%
import pandas as pd 
import glob

# Load all of the size measurements from the directory
size_data = pd.concat([pd.read_csv(f) for f in glob.glob('2024-*/processed/*sizes.csv')])

# Compute the mean cell size quantities
groups = ['date', 'replicate', 'strain', 'carbon_source', 'inducer_conc']
quants = ['width_median', 'length', 'volume', 'surface_area', 'surface_to_volume']
size_data_agged = size_data.groupby(groups)[quants].mean().reset_index()

# Rename the columns appropriately
size_data_agged.rename(columns = {'width_median':'width_um',
                                  'length':'length_um',
                                  'volume': 'volume_fL',
                                  'surface_area': 'surface_area_um2',
                                  'surface_to_volume': 'surface_to_volume_inv_um'
                                  }, inplace=True)



# Load the file that prescribes valid experiments
metadata = pd.read_csv('./valid_experimental_metadata.csv')

# Filter only on valid experiments 
valid_exps = []
for i in range(len(size_data_agged)):
    _samp = size_data_agged.iloc[i]
    

