# %%
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv(
    '../../processing/microscopy/size_measurement/output/compiled_size_measurements.csv')

data = data.groupby(['date', 'run_no', 'carbon_source', 'strain',
                     'inducer', 'inducer_conc', 'temperature_C'])[['width_median', 'length']].agg(('mean', 'sem')).reset_index()
data.to_csv(
    '../../../data/summaries/summarized_size_measurements.csv', index=False)
