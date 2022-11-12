# %%
import numpy as np
import pandas as pd
import cmdstanpy

# Load the two datasets
gr_data = pd.read_csv(
    '../../data/growth_curves/wt_growth_measurements_processed.csv')
size_data = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')
size_dates = [s[:-1] for s in size_data['date'].unique()]
gr_dates = gr_data[gr_data['date'].isin(size_dates)]

for g, d in size_data.groupby(['date', 'carbon_source']):
    _date = g[0][:-1]
    gr = gr_dates[(gr_dates['date'] == _date) & (
        gr_dates['carbon_source'] == g[1])]
    print(_date, g[1])


gr_dates[gr_dates['']]
