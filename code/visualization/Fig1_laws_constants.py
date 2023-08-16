# %%
import numpy as np
import scipy.stats
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.fluxparity
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()
const = size.fluxparity.load_constants()
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/Chure2023_collated_mass_fractions.csv')
drymass_density = pd.read_csv(
    '../../data/literature/collated_drymass_densities.csv')
protein_density = pd.read_csv(
    '../../data/literature/collated_total_protein_density.csv')
dna_prot = pd.read_csv('../../data/literature/collated_dna_protein_ratio.csv')

# %% Subset the membrane ms data
mem = ms_data[ms_data['localization'] == 'membrane'].copy()
peri = ms_data[ms_data['localization'] == 'periplasm'].copy()
mem['rho_mem'] = mem['mass_fg'] / (2 * mem['surface_area'])

# %%
# Compute the constants
mean_rho_mem = np.mean(mem['rho_mem'].values)
mean_rho_bio = np.mean(drymass_density['drymass_density_fg_fL'].values)
mean_aspect = np.mean(size_data['length_um'] / size_data['width_um'])
mean_phi_mem = np.mean(mem['mass_frac'])
mean_phi_mem_outer = np.mean(
    ms_data[ms_data['localization'] == 'membrane']['mass_frac'])
m_peri = np.mean(peri['mass_fg'])
mean_phi_dna = np.mean(dna_prot['DNA_protein_ratio'])
kappa = mean_rho_bio / mean_rho_mem
beta = 0.4558

# %%
lam_range = np.linspace(0, 2.5, 200)
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
prot_popt = scipy.stats.linregress(
    prot_data['growth_rate_hr'], np.log(prot_data['fg_protein_per_cell']))
phiRb_range = phiRb_popt[1] + phiRb_popt[0] * lam_range
prot_range = np.exp(prot_popt[1] + prot_popt[0] * lam_range)
phi_peri = m_peri / prot_range

# prefactor = 12 * mean_aspect / (3 * mean_aspect - 1)
delta = 0.024
prefactor = 12 / (3 * mean_aspect - 1)

theo = prefactor * mean_aspect * (delta + 2 * (1 + phiRb_range/beta -
                                  mean_phi_mem - phi_peri)/(kappa * mean_phi_mem))


# %%
w_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
size_data['phiRb'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width_um'] = w_popt[1] + w_popt[0] * phiRb_data['growth_rate_hr']

fig, ax = plt.subplots(1, 1, figsize=(2, 2))
for g, d in size_data.groupby('source'):
    if g != 'Si et al. 2017':
        fmt = size.viz.style_point(g)
        ax.plot(d['phiRb'], d['width_um'], **fmt)
for g, d in phiRb_data.groupby('source'):
    if g != 'Si et al. 2017':
        fmt = size.viz.style_point(g)
        ax.plot(d['mass_fraction'], d['width_um'], **fmt)
ax.plot(phiRb_range, theo, '--', color=cor['primary_red'], lw=2)
