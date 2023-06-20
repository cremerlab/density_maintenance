# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import scipy.stats
import arviz as az
mapper = size.viz.lit_mapper()
markercolors = size.viz.load_markercolors()
cor, pal = size.viz.matplotlib_style()

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ribo_data = ms_data[ms_data['localization'] == 'ribosomal sector']
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
si_chlor_data = pd.read_csv(
    '../../data/literature/Si2017/si2017_chlor_phiRb.csv')
si_clim_data = pd.read_csv(
    '../../data/literature/Si2017/si2017_clim_phiRb.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')


# %%
width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
size_data['phiRb'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width'] = width_popt[1] + \
    width_popt[0] * phiRb_data['growth_rate_hr']

# %%
mem = ms_data[ms_data['localization'] == 'membrane']
mem
rho_mem = mem['mass_fg'] / (2 * mem['surface_area'])
rho_mem


model = cmdstanpy.CmdStanModel(stan_file="constant_inference.stan")

data_dict = {'N_ms': len(mem),
             'N_size': len(size_data),
             'rho_mem': rho_mem.astype(float),
             'phi_mem': mem['mass_frac'].values.astype(float),
             'alpha': size_data['length_um'].values/size_data['width_um'].values}
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)
# %%

# %%
# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.plot(si_clim_data['phi_Rb'], si_clim_data['width_um'], mapper['Si et al., 2017']['m'], label='Si et al. 2017',
        markeredgecolor=cor['primary_black'], color=mapper['Si et al., 2017']['c'], markeredgewidth=0.25, alpha=0.75)

for g, d in ribo_data.groupby('dataset_name'):
    ax.plot(d['mass_frac'], d['width'], mapper[g]['m'], markerfacecolor=mapper[g]['c'],
            markeredgecolor=cor['primary_black'], markeredgewidth=0.25, alpha=0.75, label=g)
for g, d in size_data.groupby(['source']):
    if g == 'Si et al., 2017':
        continue
    ax.plot(d['phiRb'], d['width_um'], mapper[g]['m'], markeredgecolor=cor['primary_black'],
            markeredgewidth=0.25, color=mapper[g]['c'], alpha=0.75, label=g)
for g, d in phiRb_data.groupby(['source']):
    if g == 'Si et al. 2017':
        continue
    try:
        ax.plot(d['mass_fraction'], d['width'], markercolors[g]['m'], markeredgecolor=cor['primary_black'],
                color=markercolors[g]['c'], alpha=0.75, markeredgewidth=0.25, label=g)
    except:
        continue
ax.legend(fontsize=6, ncols=3)
# ax.
# ax[0].plot(size_data['phiRb'], size_data['width_um'], 'o')
# .plot(phiRb_data['mass_fraction'], phiRb_data['width'], 'o')
plt.ylim([0.45, 1.5])
plt.xlim([0, 0.3])
ax.set_xlabel('ribosomal mass fraction', fontsize=6)
ax.set_ylabel('average width [µm]', fontsize=6)
# %%
