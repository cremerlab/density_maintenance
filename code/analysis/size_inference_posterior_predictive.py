# %%
import matplotlib.pyplot as plt
import bokeh.io
import bebi103.viz
import corner
import bebi103
import numpy as np
import pandas as pd
import cmdstanpy
import arviz
import size.viz
cor, pal = size.viz.matplotlib_style()
cor, pal = size.viz.matplotlib_style()
_ = size.viz.bokeh_style()
bokeh.io.output_notebook()

# Load the data
data = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')
data = data[(data['width_median'] <= 2) & (data['width_median'] >= 0.2) &
            (data['length'] >= 0.2) & (data['length'] <= 7)]

# Load and compile the inferential model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/hierarchical_size_inference_uncentered.stan')

# %%
for g, d in data.groupby(['carbon_source']):
    # Assemble the data dictionary
    d['idx'] = d.groupby(['date']).ngroup() + 1
    data_dict = {'J': d['idx'].max(),
                 'N': len(d),
                 'idx': d['idx'].values.astype(int),
                 'widths': d['width_median'].values.astype(float),
                 'lengths': d['length'].values.astype(float)}

    samples = model.sample(data=data_dict)
    strain = d['strain'].values[0]
    carbon = g
    temp = d['temperature_C'].values[0]
    oe = d['overexpression'].values[0]
    ind = d['inducer_conc'].values[0]
    break
# %%
samples = arviz.from_cmdstanpy(samples)
# %%
metadata = {'strain': strain,
            'carbon': carbon,
            'temp': int(temp),
            'oe': oe,
            'ind': ind}

size.viz.diagnostic_size_viz(samples, d, metadata, './')
