# %%
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv(
    '../../processing/microscopy/size_measurement/output/compiled_size_measurements.csv')
# %%
data['delta'] = 0.0246
data.loc[data['strain'] == 'lpp14', 'delta'] = 0.0276
data.loc[data['strain'] == 'lpp21', 'delta'] = 0.0285

# %%
# Recompute dimensions
data['surface_area'] = np.pi * data['length'] * data['width_median']
data['volume'] = (np.pi / 12) * data['width_median']**2 * \
    (3 * data['length'] - data['width_median'])
data['surface_to_volume'] = data['surface_area'] / data['volume']
data['periplasm_volume'] = np.pi * data['length'] * \
    data['width_median'] * data['delta']
data['periplasm_volume_fraction'] = data['periplasm_volume'].values / \
    data['volume'].values
data['aspect_ratio'] = data['length'] / data['width_median']

# %%
# Impose bounds
data = data[(data['width_median'] >= 0.25) & (data['width_median'] <= 1.8) &
            (data['width_var'] <= 0.05)]
# %%
# Enforce consistent naming
data.loc[data['inducer'].str.lower() == 'noind', 'inducer'] = 'none'

# %%
# Drop ATC induction samples
data = data[(data['inducer'] != 'atc') & (
    data['inducer_conc'].isin([0, 2, 4, 5, 6, 10, 20, 30, 50, 100]))]

# %%
data = data.groupby(['date', 'run_no', 'carbon_source', 'strain',
                     'overexpression', 'inducer', 'inducer_conc', 'temperature_C'])[
    ['width_median', 'length', 'volume',
     'surface_area', 'surface_to_volume',
     'aspect_ratio',
     'periplasm_volume', 'periplasm_volume_fraction', 'delta']].mean().reset_index()
data.to_csv(
    '../../../data/summaries/summarized_size_measurements.csv', index=False)