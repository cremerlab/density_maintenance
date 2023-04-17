# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
import arviz as az
import cmdstanpy
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Compile the model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/literature_based_model_inference.stan')

# Load the data sets
ms_data = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')

# Aggregate the mass spec data
membrane = ms_data[ms_data['membrane'] == True]
periplasm = ms_data[ms_data['periplasm'] == True]
dfs = []
for d in [membrane, periplasm]:
    _d = d.groupby(['dataset_name', 'condition', 'growth_rate_hr'])[
        'mass_frac'].sum().reset_index()
    dfs.append(_d)
membrane, periplasm = dfs

# Assemble the data dictionary
data_dict = {
    'N_size': len(size_data),
    'N_prot': len(prot_data),
    'N_mass_spec': len(membrane),
    'delta': 0.0249,

    'widths': size_data['width_um'].values.astype(float),
    'lengths': size_data['length_um'].values.astype(float),
    'volumes': size_data['volume_um3'].values.astype(float),
    'size_lam': size_data['growth_rate_hr'].values.astype(float),

    'prot_per_cell': prot_data['fg_protein_per_cell'].values.astype(float),
    'prot_lam': prot_data['growth_rate_hr'].values.astype(float),

    'phi_mem': membrane['mass_frac'].values.astype(float),
    'phi_peri': periplasm['mass_frac'].values.astype(float),
    'ms_lam': periplasm['growth_rate_hr'].values.astype(float)
}

# %%
# Sample the model
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)
# %%

pars = ['w_min', 'w_slope', 'alpha', 'rho_prot',
        'm_peri']
fig = plt.figure(figsize=(5, 5))
fig = corner(samples, group='posterior', var_names=pars, fig=fig,
             hist_kwargs={'lw': 1}, plot_contours=False, plot_density=False, data_kwargs={'ms': 1},
             divergences=True, divergences_kwargs={'color': cor['primary_red'], 'ms': 1, 'markeredgewidth': 0})

for a in fig.axes:
    a.grid(False)


# %%
# Compute the percentiles for the parameters
upper_percs = [97.5, 87.5, 75, 62.5, 55, 50]
lower_percs = [2.5, 12.5, 25, 37.5, 45, 50]
labels = ['95%', '75%', '50%', '25%', '10%', 'median']
kwargs = {'lower_bounds': lower_percs,
          'upper_bounds': upper_percs, 'interval_labels': labels}

post = samples.posterior[pars].to_dataframe().reset_index()
post['idx'] = 1
percs = size.viz.compute_percentiles(
    post, pars, 'idx', **kwargs)
percs.to_csv(
    '../../data/mcmc/literature_model_parameter_percs.csv', index=False)

# %%
# Compute the percentiles for the size ppc
size_pars = ['w_rep', 'ell_rep', 'vol_rep']
perc_df = pd.DataFrame([])
datasets = size_data['source'].values
growth_rates = size_data['growth_rate_hr'].values
for p in size_pars:
    size_ppc = samples.posterior[p].to_dataframe().reset_index()
    size_percs = size.viz.compute_percentiles(
        size_ppc, p, f'{p}_dim_0', **kwargs)
    for i, (ds, lam) in enumerate(zip(datasets, growth_rates)):
        size_percs.loc[size_percs[f'{p}_dim_0'] == i, 'source'] = ds
        size_percs.loc[size_percs[f'{p}_dim_0'] == i, 'growth_rate_hr'] = lam
    size_percs.drop(columns=[f'{p}_dim_0'], inplace=True)
    perc_df = pd.concat([perc_df, size_percs], sort=False)
perc_df.to_csv('../../data/mcmc/literature_model_size_ppcs.csv', index=False)

# %%
# Compute the percentiles for the mass spec ppc
size_pars = ['phi_peri_rep', 'phi_mem_rep', 'rho_mem', 'rho_peri']
perc_df = pd.DataFrame([])
datasets = membrane['dataset_name'].values
growth_rates = membrane['growth_rate_hr'].values
for p in size_pars:
    ms_ppc = samples.posterior[p].to_dataframe().reset_index()
    ms_percs = size.viz.compute_percentiles(ms_ppc, p, f'{p}_dim_0', **kwargs)
    for i, (ds, lam) in enumerate(zip(datasets, growth_rates)):
        ms_percs.loc[ms_percs[f'{p}_dim_0'] == i, 'source'] = ds
        ms_percs.loc[ms_percs[f'{p}_dim_0'] == i, 'growth_rate_hr'] = lam
    ms_percs.drop(columns=[f'{p}_dim_0'], inplace=True)
    perc_df = pd.concat([perc_df, ms_percs], sort=False)
perc_df.to_csv('../../data/mcmc/literature_model_ms_ppcs.csv', index=False)

# %%
# Compute the percentiles for the protein ppcs
post = samples.posterior.prot_per_cell_rep.to_dataframe().reset_index()
percs = size.viz.compute_percentiles(
    post, 'prot_per_cell_rep', 'prot_per_cell_rep_dim_0', **kwargs)
for i, (ds, lam) in enumerate(zip(prot_data['source'].values, prot_data['growth_rate_hr'].values)):
    percs.loc[percs['prot_per_cell_rep_dim_0'] == i, 'source'] = ds
    percs.loc[percs['prot_per_cell_rep_dim_0'] == i, 'growth_rate_hr'] = lam
percs.drop(columns=['prot_per_cell_rep_dim_0'], inplace=True)
percs.to_csv(
    '../../data/mcmc/literature_model_protein_ppcs.csv', index=False)
