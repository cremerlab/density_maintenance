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
biomass = pd.read_csv(
    '../../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')
biomass
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

    'N_drymass_density': len(density),
    'drymass_density': density['drymass_density_fg_fL'].values,

    'N_drymass_total': len(biomass),
    'total_drymass': biomass['dry_mass_ug'].values,

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
plt.savefig('../../../figures/mcmc/perturbation_protein_rna_growth_ppcs.pdf')

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
plt.savefig('../../../figures/mcmc/perturbation_size_ppcs.pdf')

# %%
fig, ax = plt.subplots(2, 1, figsize=(3, 3), sharex=True)
ax[0].set_ylim([100, 400])
ax[1].set_ylim([450, 550])

# Plot the ppcs
density_ppc = samples.posterior.drymass_density_ppc.to_dataframe().reset_index()
for v in density_ppc['drymass_density_ppc'].values[::3]:
    ax[0].hlines(v, 0, 2.5, color='k', linewidth=0.1, alpha=0.25)


biomass_ppc = samples.posterior.total_drymass_ppc.to_dataframe().reset_index()
for v in biomass_ppc['total_drymass_ppc'].values[::3]:
    ax[1].hlines(v, 0, 2.5, color='k', linewidth=0.1, alpha=0.25)


for g, d in density.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)

fmt = size.viz.style_point('Basan et al. 2015')
ax[1].plot(biomass['growth_rate_hr'], biomass['dry_mass_ug'], **fmt)
ax[0].set_ylabel('drymass density [fg/fL]', fontsize=6)
ax[1].set_ylabel('total drymass [µg/OD$_{600}$mL]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.savefig('../../../figures/mcmc/perturbation_literature_ppcs.pdf')

# %%

quantities = ['prot_per_biomass',  'mem_per_biomass', 'rna_per_biomass',
              'growth_rate', 'width', 'length', 'volume',
              'surface_area', 'aspect_ratio', 'surface_to_volume',
              'phi_mem', 'rho_mem', 'phi_Rb']
units = ['ug/OD600/mL',  'ug/OD600/mL', 'ug/OD600/mL',
         'hr^-1', 'um', 'um', 'fL', 'um^2', 'dimensionless', 'um^-1',
         'dimensionless', 'fg/um^2', 'dimensionless']
post_summ = pd.DataFrame([])
mapper = {0: ['relA', 1], 1: ['meshI', 100]}
for i, q in enumerate(quantities):
    post = samples.posterior[f'{q}_mu'].to_dataframe().reset_index()
    post.dropna(inplace=True)
    for g, d in post.groupby(f'{q}_mu_dim_0'):

        mean_val = np.mean(d[f'{q}_mu'].values)
        median_val = np.median(d[f'{q}_mu'].values)
        perc = np.percentile(d[f'{q}_mu'].values, (2.5, 97.5))

        _df = pd.DataFrame({'strain': 'wildtype',
                            'overexpression': mapper[g][0],
                            'inducer_conc': mapper[g][1],
                            'carbon_source': 'glucose',
                            'quantity': q,
                            'mean_value': mean_val,
                            'median_value': median_val,
                            '2.5%': perc[0],
                            '97.5%': perc[1],
                            'unit': units[i]},
                           index=[0])
        post_summ = pd.concat([post_summ, _df], sort=False)
post_summ.to_csv('../../../data/mcmc/perturbation_posterior_parameter_summaries.csv',
                 index=False)
