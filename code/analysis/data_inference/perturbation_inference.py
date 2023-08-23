"""
Perturbation Inference
==================

"""
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./perturbation_inference.stan')

# Load the datasets
prot = pd.read_csv('../../../data/summaries/summarized_total_protein.csv')
prot = prot[prot['ug_prot_per_biomass'] > 200]
mem = pd.read_csv('../../../data/summaries/summarized_membrane_protein.csv')
rna = pd.read_csv('../../../data/summaries/summarized_total_rna.csv')
rna = rna[rna['ug_rna_per_biomass'] < 100]
sizes = pd.read_csv('../../../data/summaries/summarized_size_measurements.csv')
growth = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')
density = pd.read_csv(
    '../../../data/literature/collated_drymass_densities.csv')

mapper = {1: 1, 100: 2}
# Total protein
prot = prot[prot['overexpression'].isin(['relA', 'meshI'])].copy()
prot['idx'] = [mapper[i]
               for i in prot['inducer_conc'].values]

# Membrane protein
mem = mem[mem['overexpression'].isin(['relA', 'meshI'])]
mem = mem[mem['inducer_conc'].isin([1, 100])].copy()
mem['idx'] = [mapper[i]
              for i in mem['inducer_conc'].values]

# RNA
rna = rna[rna['overexpression'].isin(['relA', 'meshI'])]
rna = rna[rna['inducer_conc'].isin([1, 100])].copy()
rna['idx'] = [mapper[i]
              for i in rna['inducer_conc'].values]

# Sizes
sizes = sizes[sizes['overexpression'].isin(['relA', 'meshI'])]
sizes = sizes[sizes['inducer_conc'].isin([1, 100])].copy()
sizes['idx'] = [mapper[i]
                for i in sizes['inducer_conc'].values]

# Growth
growth = growth[growth['overexpression'].isin(['relA', 'meshI'])]
growth = growth[growth['inducer_conc'].isin([1, 100])].copy()
growth['idx'] = [mapper[i]
                 for i in growth['inducer_conc'].values]

# %%
data_dict = {
    'J_pert': len(mapper),

    'N_drymass': len(density),
    'drymass_density': density['drymass_density_fg_fL'].values,

    'N_sizes': len(sizes),
    'size_idx': sizes['idx'].values,
    'widths': sizes['width_median'].values,
    'lengths': sizes['length'].values,
    'surface_areas': sizes['surface_area'].values,
    'volumes': sizes['volume'].values,
    'aspect_ratios': sizes['aspect_ratio'].values,

    'N_prot': len(prot),
    'prot_idx': prot['idx'].values,
    'prot_per_biomass': prot['ug_prot_per_biomass'].values,

    'N_mem': len(mem),
    'mem_idx': mem['idx'].values,
    'mem_per_biomass': mem['ug_prot_per_biomass'].values,

    'N_rna': len(rna),
    'rna_idx': rna['idx'].values,
    'rna_per_biomass': rna['ug_rna_per_biomass'].values,

    'N_growth': len(growth),
    'growth_idx': growth['idx'].values,
    'growth_rates': growth['growth_rate_hr'].values
}

_samples = model.sample(data=data_dict, adapt_delta=0.99)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(4, 1, figsize=(2, 5), sharex=True)
ax[0].set_xlim([0, 3])
# Plot PPC distributions for non-size parameters
for i, q in enumerate(['prot', 'mem', 'rna', 'growth_rates']):
    if q != 'growth_rates':
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

titles = ['total protein', 'membrane protein']
for i, d in enumerate([prot, mem]):
    ax[i].plot(d['idx'],
               d['ug_prot_per_biomass'], '_',
               markeredgecolor=cor['primary_red'], markeredgewidth=1)
    ax[i].set_title(titles[i], fontsize=6)
    ax[i].set_ylabel('protein\n[µg / OD$_{600nm}$]', fontsize=6)

ax[2].plot(rna['idx'],
           rna['ug_rna_per_biomass'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[2].set_title('total RNA', fontsize=6)
ax[2].set_ylabel('RNA\n[µg / OD$_{600nm}$]', fontsize=6)

ax[3].plot(growth['idx'],
           growth['growth_rate_hr'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[3].set_title('growth rates', fontsize=6)
ax[3].set_ylabel('growth rate\n[hr$^{-1}$]', fontsize=6)
ax[3].set_xticks([1, 2])
ax[3].set_xticklabels(['relA, 1', 'meshI, 100'], fontsize=6)

# %%
# %% Plot PPC for size params
fig, ax = plt.subplots(5, 1, figsize=(2, 6), sharex=True)
ax[0].set_xlim([0, 3])
# Plot PPC distributions for non-size parameters
units = ['[µm]', '[µm]', '', '[µm$^2$]', '[µm$^3$]']
for i, q in enumerate(['width', 'length', 'aspect_ratio', 'surface_area', 'volume']):
    ppc = samples.posterior[f'{q}_ppc'].to_dataframe().reset_index()
    for g, d in ppc.groupby(f'{q}_ppc_dim_0'):
        ax[i].plot(1 + np.ones_like(d[::3]) * g, d[f'{q}_ppc'].values[::3] / mod,
                   '_', markeredgecolor=cor['primary_black'], markeredgewidth=0.1,
                   alpha=0.75)
    ax[i].set_title(q, fontsize=6)
    ax[i].set_ylabel(f"{q.replace('_', ' ')} {units[i]}", fontsize=6)
for g, d in sizes.groupby(['idx']):
    for i, v in enumerate(['width_median', 'length', 'aspect_ratio', 'surface_area', 'volume']):
        ax[i].plot(np.ones_like(d) * g, d[v], '_',
                   markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[4].set_xticks(list(mapper.values()))
ax[4].set_xticklabels(mapper.keys(), fontsize=6)
