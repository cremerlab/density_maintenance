# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
import cmdstanpy
import size.viz
import arviz as az
cor, pal = size.viz.matplotlib_style()

# Compile the stan model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/constancy_parameter_inference.stan')


# %%
# Load the literature datasets
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
mass_spec_data = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')
total_protein_data = pd.read_csv(
    '../../data/literature/collated_total_protein.csv')
biomass = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')

# Load our measurements as necessary and link
flow_data = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
flow_data = flow_data[flow_data['strain'] == 'wildtype']
growth_rates = pd.read_csv(
    '../../data/summaries/summarized_growth_measurements.csv')
growth_rates = growth_rates[(growth_rates['strain'] == 'wildtype') & (growth_rates['overexpression'] == 'none') &
                            (growth_rates['inducer_conc'] == 0)]
growth_rates = growth_rates.groupby(['carbon_source'])[
    'growth_rate_hr'].mean().reset_index()
growth_rates_dict = {k: v for k, v in zip(growth_rates['carbon_source'].values,
                                          growth_rates['growth_rate_hr'].values)}
for k, v in growth_rates_dict.items():
    flow_data.loc[flow_data['carbon_source'] == k, 'growth_rate_hr'] = v
flow_data


# Map components of mass spec data
mass_spec_data.loc[mass_spec_data['go_terms'].str.contains(
    'GO:0005886'), 'classification'] = 'membrane'
mass_spec_data.loc[mass_spec_data['go_terms'].str.contains(
    'GO:0042597'), 'classification'] = 'periplasm'
mass_spec_data = mass_spec_data[~mass_spec_data['classification'].isnull()]
mass_spec_data = mass_spec_data.groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'classification'])['mass_frac'].sum().reset_index()

# Drop a bad mori point
mass_spec_data = mass_spec_data[~((mass_spec_data['dataset_name'] == 'Mori et al. 2021') & (
    mass_spec_data['condition'] == '0.2% glucose') & (mass_spec_data['growth_rate_hr'].isin([0.56, 0.69])))]
membrane = mass_spec_data[mass_spec_data['classification'] == 'membrane']
periplasm = mass_spec_data[mass_spec_data['classification'] == 'periplasm']


# %%
# Assemble the data dictionary
data_dict = {
    'N_mass_spec': len(membrane),
    'N_size': len(lit_size_data),
    'N_flow': len(flow_data),
    'N_prot': len(total_protein_data),
    'N_biomass': len(biomass),

    'phi_peri': periplasm['mass_frac'].values.astype(float),
    'phi_memb': membrane['mass_frac'].values.astype(float),
    'mass_spec_lambda': periplasm['growth_rate_hr'].values.astype(float),

    'size_lambda': lit_size_data['growth_rate_hr'].values.astype(float),
    'width': lit_size_data['width_um'].values.astype(float),
    'length': lit_size_data['length_um'].values.astype(float),

    'flow_meas': flow_data['cells_per_biomass'].values.astype(float),
    'flow_lambda': flow_data['growth_rate_hr'].values.astype(float),

    'prot_meas': total_protein_data['total_protein_ug_od600'].values.astype(float),
    'prot_lambda': total_protein_data['growth_rate_hr'].values.astype(float),
    'biomass': biomass['dry_mass_ug'].values.astype(float)
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)
# %%
model_params = samples.posterior[[
    'm_peri', 'rho_mem', 'alpha']].to_dataframe().reset_index()
model_params = model_params[['m_peri', 'rho_mem', 'alpha']]
fig = plt.figure(figsize=(3, 3))
fig = corner(model_params, color=cor['primary_blue'], fig=fig,
             hist_kwargs={'lw': 1}, plot_contours=False, plot_density=False, data_kwargs={'ms': 1, 'color': cor['primary_blue']})

for a in fig.axes:
    a.grid(False)
# %%
mapper = size.viz.lit_mapper()
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
ax[0].set_ylim([0.01, 0.12])
ax[1].set_ylim([0.01, 0.20])
ax[2].set_ylim([0, 1])
for g in mass_spec_data['dataset_name'].unique():
    p = periplasm[periplasm['dataset_name'] == g]
    m = membrane[membrane['dataset_name'] == g]
    ax[0].plot(p['growth_rate_hr'].values, p['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
    ax[1].plot(m['growth_rate_hr'].values, m['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
    ax[2].plot(p['growth_rate_hr'].values, p['mass_frac'].values/m['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
