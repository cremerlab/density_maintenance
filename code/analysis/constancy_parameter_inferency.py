# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
import cmdstanpy
import size.viz
import arviz as az
import imp
imp.reload(size.viz)
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# Compile the stan model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/growth_rate_dependence_inference.stan')

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
mass_spec_data = mass_spec_data[~((mass_spec_data['dataset_name'] == 'Mori et al. 2021') & (
    mass_spec_data['condition'] == '0.2% glucose') & (mass_spec_data['growth_rate_hr'].isin([0.56, 0.69])))]
membrane = mass_spec_data[(mass_spec_data['classification'] == 'membrane')]  # & (
# mass_spec_data['periplasm'] == False)]
periplasm = mass_spec_data[(mass_spec_data['classification'] == 'periplasm')]  # & (
# mass_spec_data['periplasm'] == True)]
membrane = membrane.groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'classification'])['mass_frac'].sum().reset_index()
periplasm = periplasm.groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'classification'])['mass_frac'].sum().reset_index()


# %%
# Assemble the data dictionary
data_dict = {
    'const_phi_mem': 1,
    'N_mass_spec': len(membrane),
    'N_size': len(lit_size_data),
    'N_flow': len(flow_data),
    'N_prot': len(total_protein_data),
    'N_biomass': len(biomass),

    'phi_peri': periplasm['mass_frac'].values.astype(float),
    'phi_mem': membrane['mass_frac'].values.astype(float),
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
# model_params = samples.posterior[[
# 'm_peri', 'rho_mem', 'alpha', 'kappa']].to_dataframe().reset_index()
# pars = ['m_peri', 'rho_mem', 'alpha', 'kappa', 'k_m']

pars = ['w_min', 'alpha', 'm_min', 'k_w',
        'k_m', 'm_peri', 'phi_mem_mu']
fig = plt.figure(figsize=(5, 5))
fig = corner(samples, group='posterior', var_names=pars, fig=fig,
             hist_kwargs={'lw': 1}, plot_contours=False, plot_density=False, data_kwargs={'ms': 1},
             divergences=True, divergences_kwargs={'color': cor['primary_red'], 'ms': 1, 'markeredgewidth': 0})

for a in fig.axes:
    a.grid(False)


# %%
# lam_range = np.linspace(0.05, 3, 100)
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 7])
ax[2].set_ylim([0, 5])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]')
ax[0].set_ylabel('width [µm]')
ax[1].set_ylabel('length [µm]')
ax[2].set_ylabel('volume [µm]')
ax[3].set_ylabel('protein / biomass\n[µg / ODmL]')
# ax[2].set_xlim([0, 1.75])
# ax[2].set_ylim([275, 450])

pars = ['width_rep', 'length_rep', 'volume_rep', 'prot_rep']
for j, p in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    perc = size.viz.compute_percentiles(post, p, f'{p}_dim_0')
    if j <= 2:
        d = lit_size_data
    else:
        d = total_protein_data
    for i in range(len(d)):
        perc.loc[perc[f'{p}_dim_0'] == i,
                 'source'] = d['source'].values[i]
        perc.loc[perc[f'{p}_dim_0'] == i,
                 'growth_rate_hr'] = d['growth_rate_hr'].values[i]

    for g, d in perc.groupby(['interval']):
        d.sort_values(by='growth_rate_hr', inplace=True)
        ax[j].fill_between(d['growth_rate_hr'], d['lower'],
                           d['upper'], color=cor['blue'], alpha=0.1)

    for g, d in lit_size_data.groupby(['source']):
        ax[0].plot(d['growth_rate_hr'], d['width_um'], mapper[g]['m'], ms=5, markerfacecolor=mapper[g]['c'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, alpha=0.75)
        ax[1].plot(d['growth_rate_hr'], d['length_um'], mapper[g]['m'], ms=5, markerfacecolor=mapper[g]['c'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, alpha=0.75)
        ax[2].plot(d['growth_rate_hr'], d['volume_um3'], mapper[g]['m'], ms=5, markerfacecolor=mapper[g]['c'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, alpha=0.75)

for g, d in total_protein_data.groupby('source'):
    ax[3].plot(d['growth_rate_hr'], d['total_protein_ug_od600'], mapper[g]['m'], ms=5, markerfacecolor=mapper[g]['c'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.5, alpha=0.75)

plt.tight_layout()

# %%
fig, ax = plt.subplots(1, 3, figsize=(8, 2))
ax[0].set_ylim([0, 0.25])
# ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 1])
# ax[1].set_ylim([0, 15])

vars = ['phi_mem_rep', 'phi_peri_rep', 'rel_phi_rep']
for i, v in enumerate(vars):
    _df = samples.posterior[f'{v}'].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(_df, v, f'{v}_dim_0')
    for j, k in enumerate(membrane['growth_rate_hr'].values):
        percs.loc[percs[f'{v}_dim_0'] == j, 'growth_rate_hr'] = k
    for g, d in percs.groupby(['interval'], sort=False):
        d.sort_values(by='growth_rate_hr', inplace=True)
        ax[i].fill_between(d['growth_rate_hr'], d['lower'],
                           d['upper'], alpha=0.2, color=cor['blue'])

# kappa_df = samples.posterior.kappa.to_dataframe().reset_index()
# kappa_df['idx'] = 1
# percs = size.viz.compute_percentiles(kappa_df, 'kappa', 'idx')


for g in mass_spec_data['dataset_name'].unique():
    p = periplasm[periplasm['dataset_name'] == g]
    m = membrane[membrane['dataset_name'] == g]
    ax[0].plot(p['growth_rate_hr'].values, m['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
    ax[1].plot(p['growth_rate_hr'].values, p['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
    ax[2].plot(p['growth_rate_hr'].values, p['mass_frac'].values/m['mass_frac'].values, linestyle='none', marker=mapper[g]['m'],
               color=mapper[g]['c'], alpha=0.5, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)

for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]')
ax[0].set_ylabel('membrane protein mass fraction')
ax[1].set_ylabel('periplasmic protein mass fraction')
ax[2].set_ylabel('relative mass fraction')
