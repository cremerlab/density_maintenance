"""
Wildtype Inference
==================

This script executes a bayesian inferential model to estimate size and protein 
load parameters from measurements of wildtype E. coli.  This generates figures 
illustrating the posterior predictive checks of the size, growth, and protein 
inference as well as CSV files of the summarized posteriors for each parameter. 
"""
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./wildtype_inference.stan')

# Load the datasets
# prot = pd.read_csv('../../../data/summaries/summarized_total_protein.csv')
# rna = pd.read_csv('../../../data/summaries/summarized_total_rna.csv')
phi = pd.read_csv(
    '../../../data/summaries/merged_rna_protein_measurements.csv')
mem = pd.read_csv('../../../data/summaries/summarized_membrane_protein.csv')
peri = pd.read_csv(
    '../../../data/summaries/summarized_periplasmic_protein.csv')
peri = peri[peri['ug_prot_per_biomass'] <= 50]
sizes = pd.read_csv('../../../data/summaries/summarized_size_measurements.csv')
flow = pd.read_csv('../../../data/summaries/summarized_cell_counts.csv')
flow = flow[flow['cells_per_biomass'].values > 1E6]
flow = flow[flow['strain'] == 'wildtype']
growth = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')

# %%
# Restrict to wildtype
phi = phi[(phi['strain'] == 'wildtype') &
          (phi['overexpression'] == 'none')]
mem = mem[(mem['strain'] == 'wildtype') & (mem['overexpression'] == 'none')]
peri = peri[(peri['strain'] == 'wildtype') &
            (peri['overexpression'] == 'none')]
sizes = sizes[(sizes['strain'] == 'wildtype') &
              (sizes['overexpression'] == 'none')]
sizes = sizes[sizes['carbon_source'] != 'ezMOPS']
growth = growth[(growth['strain'] == 'wildtype') &
                (growth['overexpression'] == 'none')]

# Add a condition mapper
mapper = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3,
          'glucose': 4, 'glucoseCAA': 5, 'LB': 6}

for d in [phi, mem, peri, sizes, growth, flow]:
    d['idx'] = [mapper[k] for k in d['carbon_source'].values]
    d['idx'] = d['idx'].values.astype(int)

data_dict = {
    'N_phi': len(phi),
    'J_phi': phi['idx'].max(),
    'phi_idx': phi['idx'].values,
    'prot_per_biomass': phi['ug_prot_per_biomass'].values,
    'rna_per_biomass': phi['ug_rna_per_biomass'].values,
    'phiRb': phi['phiRb'].values,

    'N_peri': len(peri),
    'J_peri': peri['idx'].max(),
    'peri_idx': peri['idx'].values,
    'peri_per_biomass': peri['ug_prot_per_biomass'].values,

    'N_mem': len(mem),
    'J_mem': mem['idx'].max(),
    'mem_idx': mem['idx'].values,
    'mem_per_biomass': mem['ug_prot_per_biomass'].values,

    'N_growth': len(growth),
    'J_growth': growth['idx'].max(),
    'growth_idx': growth['idx'].values,
    'growth_rates': growth['growth_rate_hr'].values,

    'N_size': len(sizes),
    'J_size': sizes['idx'].max(),
    'size_idx': sizes['idx'].values,
    'widths': sizes['width_median'].values,
    'lengths': sizes['length'].values,
    'volumes': sizes['volume'].values,
    'surface_areas': sizes['surface_area'].values,
    'surface_to_volumes': sizes['surface_to_volume'].values,
    'aspect_ratios': sizes['aspect_ratio'].values,

    'N_flow': len(flow),
    'J_flow': flow['idx'].max(),
    'flow_idx': flow['idx'].values,
    'cells_per_biomass': flow['cells_per_biomass'].values,
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(7, 1, figsize=(3, 8), sharex=True)
# Plot PPC distributions for non-size parameters
for i, q in enumerate(['prot', 'mem', 'peri', 'rna', 'phiRb', 'cells', 'growth_rates']):
    if (q != 'growth_rates') & (q != 'phiRb'):
        n = f'{q}_per_biomass'
    elif q == 'phiRb':
        n = q
    else:
        n = 'growth_rates'
    if q == 'cells':
        mod = 1E9
    else:
        mod = 1
    ppc = samples.posterior[f'{n}_ppc'].to_dataframe().reset_index()

    for g, d in ppc.groupby(f'{n}_ppc_dim_0'):
        ax[i].plot(1 + np.ones_like(d[::3]) * g, d[f'{n}_ppc'].values[::3] / mod,
                   '_', markeredgecolor=cor['primary_black'], markeredgewidth=0.1,
                   alpha=0.75)

titles = ['total protein', 'membrane protein',  'periplasmic protein']
for i, d in enumerate([phi, mem, peri]):
    ax[i].plot(d['idx'],
               d['ug_prot_per_biomass'], '_',
               markeredgecolor=cor['primary_red'], markeredgewidth=1)
    ax[i].set_title(titles[i], fontsize=6)
    ax[i].set_ylabel('protein\n[µg / OD$_{600nm}$]', fontsize=6)

ax[3].plot(phi['idx'],
           phi['ug_rna_per_biomass'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[3].set_title('total RNA', fontsize=6)
ax[3].set_ylabel('RNA\n[µg / OD$_{600nm}$]', fontsize=6)

ax[4].plot(phi['idx'],
           phi['phiRb'], '_',
           markeredgecolor=cor['primary_red'],
           markeredgewidth=0.5)
ax[4].set_title('ribosomal allocation', fontsize=6)
ax[4].set_ylabel('$\phi_{Rb}$', fontsize=6)

ax[5].plot(flow['idx'],
           flow['cells_per_biomass'] / 1E9, '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[5].set_title('flow cytometry', fontsize=6)
ax[5].set_ylabel('cells per biomass\n' + r'[$\times 10^9$]', fontsize=6)

ax[6].plot(growth['idx'],
           growth['growth_rate_hr'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[6].set_title('growth rates', fontsize=6)
ax[6].set_ylabel('growth rate\n[hr$^{-1}$]', fontsize=6)
ax[6].set_xticks(list(mapper.values()))
ax[6].set_xticklabels(mapper.keys(), fontsize=6)
plt.savefig('../../../figures/mcmc/wildtype_protein_flow_growth_ppcs.pdf',
            bbox_inches='tight')

# %% Plot PPC for size params
fig, ax = plt.subplots(6, 1, figsize=(3, 6), sharex=True)
# Plot PPC distributions for non-size parameters
units = ['[µm]', '[µm]', '', '[µm$^2$]', '[µm$^3$]', '[µm$^{-1}$]']
for i, q in enumerate(['width', 'length', 'aspect_ratio', 'surface_area', 'volume', 'surface_to_volume']):
    ppc = samples.posterior[f'{q}_ppc'].to_dataframe().reset_index()
    for g, d in ppc.groupby(f'{q}_ppc_dim_0'):
        ax[i].plot(1 + np.ones_like(d[::3]) * g, d[f'{q}_ppc'].values[::3] / mod,
                   '_', markeredgecolor=cor['primary_black'], markeredgewidth=0.1,
                   alpha=0.75)
    ax[i].set_title(q, fontsize=6)
    ax[i].set_ylabel(f"{q.replace('_', ' ')} {units[i]}", fontsize=6)
for g, d in sizes.groupby(['carbon_source']):
    for i, v in enumerate(['width_median', 'length', 'aspect_ratio', 'surface_area', 'volume', 'surface_to_volume']):
        ax[i].plot(np.ones_like(d) * mapper[g], d[v], '_',
                   markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[5].set_xticks(list(mapper.values()))
ax[5].set_xticklabels(mapper.keys(), fontsize=6)
plt.savefig('../../../figures/mcmc/wildtype_size_ppcs.pdf',
            bbox_inches='tight')

# %%
# Summarize the posteriors with mean and 95th percentile
quantities = ['prot_per_biomass', 'peri_per_biomass', 'prot_per_cell', 'mem_per_biomass', 'rna_per_biomass',
              'cells_per_biomass', 'growth_rate', 'width', 'length', 'volume',
              'surface_area', 'aspect_ratio', 'surface_to_volume', 'm_peri',
              'rho_peri', 'phi_peri', 'phi_mem', 'rho_mem', 'phiRb']
units = ['ug/OD600/mL', 'ug/OD600/mL', 'fg/cell', 'ug/OD600/mL', 'ug/OD600/mL',
         'cells/OD600/mL',  'hr^-1', 'um', 'um', 'fL', 'um^2', 'dimensionless', 'um^-1',
         'fg/cell', 'fg/fL', 'dimensionless', 'dimensionless', 'fg/um^2', 'dimensionless']
post_summ = pd.DataFrame([])
rev_mapper = {v: k for k, v in mapper.items()}
for i, q in enumerate(quantities):
    post = samples.posterior[f'{q}_mu'].to_dataframe().reset_index()
    post.dropna(inplace=True)
    for g, d in post.groupby(f'{q}_mu_dim_0'):

        mean_val = np.mean(d[f'{q}_mu'].values)
        median_val = np.median(d[f'{q}_mu'].values)
        perc = np.percentile(d[f'{q}_mu'].values, (2.5, 97.5))

        _df = pd.DataFrame({'strain': 'wildtype',
                            'overexpression': 'none',
                            'inducer_conc': 0,
                            'carbon_source': rev_mapper[g+1],
                            'quantity': q,
                            'mean_value': mean_val,
                            'median_value': median_val,
                            '2.5%': perc[0],
                            '97.5%': perc[1],
                            'unit': units[i]},
                           index=[0])
        post_summ = pd.concat([post_summ, _df], sort=False)
post_summ.to_csv('../../../data/mcmc/wildtype_posterior_parameter_summaries.csv',
                 index=False)
