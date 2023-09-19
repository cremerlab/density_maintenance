# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import size.viz
import matplotlib.pyplot as plt
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./lpp_inference.stan')

size_data = pd.read_csv(
    '../../../data/summaries/summarized_size_measurements.csv')
size_data = size_data[size_data['strain'] == 'lpp14']
size_data = size_data[size_data['carbon_source'] != 'RDM']
flow_data = pd.read_csv('../../../data/summaries/summarized_cell_counts.csv')
flow_data = flow_data[flow_data['strain'] == 'lpp14']
peri_data = pd.read_csv(
    '../../../data/summaries/summarized_periplasmic_protein.csv')
peri_data = peri_data[peri_data['strain'] == 'lpp14']
growth_data = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')
growth_data = growth_data[growth_data['strain'] == 'lpp14']
growth_data = growth_data[growth_data['carbon_source'] != 'ezMOPS']
# %%
# Add mapper
mapper = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3,
          'glucose': 4, 'glucoseCAA': 5, 'LB': 6}
for d in [size_data, flow_data, peri_data, growth_data]:
    for k, v in mapper.items():
        d.loc[d['carbon_source'] == k, 'idx'] = v
    d['idx'] = d['idx'].values.astype(int)

# %%
# Set up the data dictionary
data_dict = {'N_peri': len(peri_data),
             'J_peri': peri_data['idx'].max(),
             'peri_idx': peri_data['idx'].values,
             'peri_per_biomass': peri_data['ug_prot_per_biomass'].values,

             'N_growth': len(growth_data),
             'J_growth': growth_data['idx'].max(),
             'growth_idx': growth_data['idx'].values,
             'growth_rates': growth_data['growth_rate_hr'].values,

             'N_size': len(size_data),
             'J_size': size_data['idx'].max(),
             'size_idx': size_data['idx'].values,
             'widths': size_data['width_median'].values,
             'lengths': size_data['length'].values,
             'aspect_ratios': size_data['aspect_ratio'].values,
             'surface_areas': size_data['surface_area'].values,
             'surface_to_volumes': size_data['surface_to_volume'].values,
             'volumes': size_data['volume'].values,

             'N_flow': len(flow_data),
             'J_flow': flow_data['idx'].max(),
             'flow_idx': flow_data['idx'].values,
             'cells_per_biomass': flow_data['cells_per_biomass'].values}

_samples = model.sample(data_dict, adapt_delta=0.9999)
samples = az.from_cmdstanpy(_samples)

# %%o

fig, ax = plt.subplots(3, 1, figsize=(3, 6), sharex=True)
# Plot PPC distributions for non-size parameters
for i, q in enumerate(['peri', 'cells', 'growth_rates']):
    if (q != 'growth_rates') & (q != 'phiRb'):
        n = f'{q}_per_biomass'
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
ax[0].plot(peri_data['idx'], peri_data['ug_prot_per_biomass'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[0].set_title(titles[i], fontsize=6)
ax[0].set_ylabel('protein\n[µg / OD$_{600nm}$]', fontsize=6)

ax[1].plot(flow_data['idx'],
           flow_data['cells_per_biomass'] / 1E9, '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[1].set_title('flow cytometry', fontsize=6)
ax[1].set_ylabel('cells per biomass\n' + r'[$\times 10^9$]', fontsize=6)

ax[2].plot(growth_data['idx'],
           growth_data['growth_rate_hr'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[2].set_title('growth rates', fontsize=6)
ax[2].set_ylabel('growth rate\n[hr$^{-1}$]', fontsize=6)
ax[2].set_xticks(list(mapper.values()))
ax[2].set_xticklabels(mapper.keys(), fontsize=6)
plt.savefig('../../../figures/mcmc/lpp_protein_flow_growth_ppcs.pdf',
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
for g, d in size_data.groupby(['carbon_source']):
    for i, v in enumerate(['width_median', 'length', 'aspect_ratio', 'surface_area', 'volume', 'surface_to_volume']):
        ax[i].plot(np.ones_like(d) * mapper[g], d[v], '_',
                   markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[5].set_xticks(list(mapper.values()))
ax[5].set_xticklabels(mapper.keys(), fontsize=6)
plt.savefig('../../../figures/mcmc/wildtype_size_ppcs.pdf',
            bbox_inches='tight')
