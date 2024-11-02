#%%
import pandas as pd 
import glob
import tqdm

#%% Load the mass spectrometry data
data = pd.read_csv('./mass_spectrometry_allocation_wide.csv')
data = data.select_dtypes(include=['float64']).round(decimals=3).combine_first(data)

groupby = ['date', 'replicate', 'strain', 'carbon_source', 'temperature_C',
           'inducer_conc']
merged = pd.DataFrame([])
for date in tqdm.tqdm(data['date'].unique()):
    # Parse the valid growth and size files
    growth_files = glob.glob(f'./{date}_r*/processed/*growth_rates.csv') 
    size_files = glob.glob(f'./{date}_r*/processed/*sizes.csv')
    for g, s in zip(growth_files, size_files):
        growth = pd.read_csv(g)
        sizes = pd.read_csv(s)

        # Handle edge case for underspecified data. 
        for _d in [growth, sizes]:
            if 'strain' not in _d.keys():
                _d['strain'] = 'wildtype'
            if 'temperature_C' not in _d.keys():
                _d['temperature_C'] = 37
            if 'inducer_conc' not in _d.keys():
                _d['inducer_conc'] = 0

        # Aggregate size data
        sizes = sizes.groupby(groupby)[['length', 'width_median', 'volume', 
                                        'surface_area', 'surface_to_volume',
                                        'periplasm_volume']].mean().reset_index()  

        # Rename the size cols
        sizes.rename(columns={'length':'length_um', 'width_median':'width_um',
                'volume':'volume_fL', 'surface_area':'surface_area_um2',
                'surface_to_volume':'sav_inv_um'}, inplace=True)

        # Round sizes and growth rates to three decimal places
        sizes = sizes.select_dtypes(include=['float64']
                                    ).round(decimals=3).combine_first(sizes)
        growth = growth.select_dtypes(include=['float64']).round(decimals=3).combine_first(growth)
        growth.drop(columns=['od_init', 'growth_rate_std'], inplace=True)
        growth['temperature_C'] = 37 # Add in missing data. 

        # Merge the two
        _merged = growth.merge(sizes, on=groupby, how='outer')
        merged = pd.concat([merged, _merged]) 
# Drop any continuance of lacZ data at this point
merged = merged[merged['strain'] != 'lacZ']
#%%
data['temperature_C'] = 37
joined = data.merge(merged, on=groupby, how='outer')

# Drop one invalid case (pRelA + 2ng/ml dox, rep 2 from 2024-05-23), the only 
# nan
joined.dropna(inplace=True)

# Save as complete dataset 
joined.to_csv('../../../data/compiled_measurements.csv', index=False)