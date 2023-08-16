# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import arviz as az
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# Load the necessary data sets
# growth_data = pd.read_csv('../../data/processed/labeled_growth_measurements.csv')
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
lit_size_data = lit_size_data[(lit_size_data['source']
                              != 'Zaritsky et al. 1993') &
                              (lit_size_data['source'] != 'Grossman et al. 1982') &
                              (lit_size_data['source'] != 'Zaritsky & Woldringh 1978')]

model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/hierarchical_envelope_growth_laws.stan')

# %%
lit_size_data['idx'] = lit_size_data.groupby('source').ngroup() + 1
data_dict = {'N_size': len(lit_size_data),
             'J_size': lit_size_data['idx'].max(),
             'size_idx': lit_size_data['idx'].values,
             'size_lam': lit_size_data['growth_rate_hr'].values,
             'aspect_ratio': lit_size_data['length_um'].values / lit_size_data['width_um'].values}
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(lit_size_data['idx'].max(
), 1, figsize=(4, 10), sharey=True, sharex=True)
for g, d in lit_size_data.groupby(['source', 'idx']):
    fmt = size.viz.style_point(g[0])
    ax[g[1] - 1].plot(d['growth_rate_hr'],
                      d['length_um'] / d['width_um'], **fmt)

ax[-1].legend()
# %%
df = samples.posterior.alpha_mu.to_dataframe().reset_index()
plt.hist(df['alpha_mu'], bins=20)
