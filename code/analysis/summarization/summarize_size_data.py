# %%
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv(
    '../../processing/microscopy/size_measurement/output/compiled_size_measurements.csv')
data['periplasm_volume_fraction'] = data['periplasm_volume'].values / data['volume'].values
# Impose bounds
data = data[(data['width_median'] >= 0.3) & (data['width_median'] <= 1.5) &
            (data['surface_to_volume'] <= 8)]
# %%
data = data.groupby(['date', 'run_no', 'carbon_source', 'strain',
                     'inducer', 'inducer_conc', 'temperature_C'])[
    ['width_median', 'length', 'volume',
     'surface_area', 'surface_to_volume', 'periplasm_volume', 'periplasm_volume_fraction']].mean().reset_index()
data.to_csv(
    '../../../data/summaries/summarized_size_measurements.csv', index=False)

#%%