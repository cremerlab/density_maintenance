# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()

size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']

prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
mem_data = ms_data[ms_data['localization'] == 'membrane']
peri_data = ms_data[ms_data['localization'] == 'periplasm']
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']

# %%
model = cmdstanpy.CmdStanModel(stan_file='mass_spec_empirics_uncertainty.stan')

# %%
data_dict = {
    'N_size': len(size_data),
    'size_lam': size_data['growth_rate_hr'].values,
    'widths': size_data['width_um'].values,
    'lengths': size_data['length_um'].values,

    'N_prot': len(prot_data),
    'prot_lam': prot_data['growth_rate_hr'].values,
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,

    'N_ms': len(mem_data),
    'ms_lam': mem_data['growth_rate_hr'].values,
    'phi_mem': mem_data['mass_frac'].values,
    'phi_peri': peri_data['mass_frac'].values,

    'N_phi': len(phiRb_data),
    'phi_lam': phiRb_data['growth_rate_hr'].values
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
# Plot ppcs
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_xlim([0.01, 2.5])
ax[0].set_ylabel('width [µm]', fontsize=6)
ax[1].set_ylabel('length [µm]', fontsize=6)
ax[2].set_ylabel('protein [fg / cell]', fontsize=6)

# PPCs
width_ppcs = samples.posterior.width_ppc.to_dataframe().reset_index()
length_ppcs = samples.posterior.length_ppc.to_dataframe().reset_index()
prot_ppcs = samples.posterior.prot_per_cell_ppc.to_dataframe().reset_index()

for i, (g, d) in enumerate(width_ppcs.groupby(['chain', 'draw'])):
    if i % 10 == 0:
        ax[0].plot(size_data['growth_rate_hr'],
                   d['width_ppc'], 'k-', lw=0.1, alpha=0.1)

for i, (g, d) in enumerate(length_ppcs.groupby(['chain', 'draw'])):
    if i % 10 == 0:
        ax[1].plot(size_data['growth_rate_hr'],
                   d['length_ppc'], 'k-', lw=0.1, alpha=0.1)

for i, (g, d) in enumerate(prot_ppcs.groupby(['chain', 'draw'])):
    if i % 10 == 0:
        ax[2].plot(prot_data['growth_rate_hr'],
                   d['prot_per_cell_ppc'], 'k-', lw=0.1, alpha=0.1)


# Data
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['width_um'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], **fmt)

for g, d in prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

# %%
sources = [m for m in mem_data['dataset_name'].values]
lam = [ell for ell in mem_data['growth_rate_hr'].values]

fig, ax = plt.subplots(2, 1, figsize=(4, 4), sharex=True)

# Generate the ppcs
peri_ppcs = samples.posterior.m_peri.to_dataframe().reset_index()
rho_ppcs = samples.posterior.rho_mem.to_dataframe().reset_index()
for g, d in peri_ppcs.groupby(['m_peri_dim_0']):
    lower, upper = np.percentile(d['m_peri'].values, (2.5, 97.5))
    median = np.median(d['m_peri'].values)
    fmt = size.viz.style_point(sources[g])
    ax[0].vlines(lam[g], lower, upper, lw=0.75, alpha=0.5, color=fmt['color'])
    ax[0].plot(lam[g], median, **fmt)

for g, d in rho_ppcs.groupby(['rho_mem_dim_0']):
    lower, upper = np.percentile(d['rho_mem'].values, (2.5, 97.5))
    median = np.median(d['rho_mem'].values)
    fmt = size.viz.style_point(sources[g])
    ax[1].vlines(lam[g], lower, upper, lw=1, alpha=0.5, color=fmt['color'])
    ax[1].plot(lam[g], median, **fmt)
ax[0].set_ylim([0, 20])
ax[1].set_ylim([0, 6])

ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('periplasmic protein [fg / cell]', fontsize=6)
ax[1].set_ylabel('membrane protein density [fg / µm$^2$]', fontsize=6)

# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

# for g, d in phiRb_data.groupby('source'):
#     fmt = size.viz.style_point(g)
#     ax.plot(d['growth_rate_hr'], d['mass_fraction'], **fmt)

sources = [s for s in phiRb_data['source'].values]
phis = [phi for phi in phiRb_data['mass_fraction'].values]
ppcs = samples.posterior.phi_width.to_dataframe().reset_index()

for g, d in ppcs.groupby(['phi_width_dim_0']):
    fmt = size.viz.style_point(sources[g])
    lower, upper = np.percentile(d['phi_width'].values, (2.5, 97.5))
    median = np.median(d['phi_width'].values)
    ax.vlines(phis[g], lower, upper, lw=0.75, alpha=0.5, color=fmt['color'])
    ax.plot(phis[g], median, **fmt)

ax.set_xlabel('allocation towards ribosomes', fontsize=6)
ax.set_ylabel('cell width [µm]', fontsize=6)
