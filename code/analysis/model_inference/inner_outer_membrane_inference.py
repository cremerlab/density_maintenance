# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import size.viz

size_data = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
inner_mem = ms_data[ms_data['localization'] == 'inner membrane']
outer_mem = ms_data[ms_data['localization'] == 'outer membrane']
prot_data = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')

# Load the inferential model
model = cmdstanpy.CmdStanModel(stan_file='inner_outer_membrane_inference.stan')

# Define the data dictionary
data_dict = {
    'N_size': len(size_data),
    'size_lam': size_data['growth_rate_hr'].values,
    'surface_areas': size_data['surface_area_um2'].values,

    'N_prot': len(prot_data),
    'prot_lam': prot_data['growth_rate_hr'].values,
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,

    'N_ms': len(inner_mem),
    'ms_lam': inner_mem['growth_rate_hr'].values,
    'phi_mem_outer': outer_mem['mass_frac'].values,
    'phi_mem_inner': inner_mem['mass_frac'].values,
}

_samples = model.sample(data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
pars = ['ms_rho_mem_inner', 'ms_rho_mem_outer']
percs = pd.DataFrame([])
lower = [2.5, 12.5, 37.5, 50]
upper = [97.5, 87.5, 62.5, 50]
labels = ['95%', '75%', '25%', 'median']
for p in pars:
    post = samples.posterior[p].to_dataframe().reset_index()
    perc_df = size.viz.compute_percentiles(
        post, p, groupby=f'{p}_dim_0', lower_bounds=lower, upper_bounds=upper, interval_labels=labels)
    perc_df['source'] = [inner_mem['dataset_name'].values[k]
                         for k in perc_df[f'{p}_dim_0']]
    perc_df['growth_rate_hr'] = [inner_mem['growth_rate_hr'].values[k]
                                 for k in perc_df[f'{p}_dim_0']]
    perc_df = perc_df[['source', 'growth_rate_hr',
                       'quantity', 'interval', 'lower', 'upper']]
    percs = pd.concat([percs, perc_df], sort=False)
percs.to_csv(
    '../../../data/mcmc/inner_outer_membrane_densities_longform.csv', index=False)

# %%
wide_percs = pd.DataFrame([])
for g, d in percs.groupby(['quantity', 'source', 'growth_rate_hr']):
    _df = pd.DataFrame({'source': g[1],
                        'quantity': g[0],
                        'growth_rate_hr': g[2],
                        'median_value': d[d['interval'] == 'median']['lower'].values[0],
                        '2.5%': d[d['interval'] == '95%']['lower'].values[0],
                        '97.5%': d[d['interval'] == '95%']['upper'].values[0]},
                       index=[0])
    wide_percs = pd.concat([wide_percs, _df], sort=False)
wide_percs.to_csv(
    '../../../data/mcmc/inner_outer_membrane_densities_wide.csv', index=False)
