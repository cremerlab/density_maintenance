# %%
import corner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import size.fluxparity
import scipy.stats
import arviz as az
mapper = size.viz.lit_mapper()
markercolors = size.viz.load_markercolors()
consts = size.fluxparity.load_constants()
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
rho_mem = mem['mass_fg'] / (2 * mem['surface_area'])
N_ppc = 500
nu_range = np.linspace(0.1, 20, N_ppc)
phiRb_range = size.fluxparity.phiRb_optimal_allocation(consts['gamma_max'], nu_range,
                                                       consts['Kd_cpc'], consts['phi_O'])
lam_range = size.fluxparity.steady_state_growth_rate(
    consts['gamma_max'], phiRb_range, nu_range, consts['Kd_cpc'], consts['phi_O'])

model = cmdstanpy.CmdStanModel(stan_file="constant_inference.stan")

data_dict = {'N_ms': len(mem),
             'N_size': len(size_data),
             'N_ppc': N_ppc,
             'phiRb': phiRb_range,
             'delta': 0.0246,
             'beta_rp': 0.4558,
             'rho_mem': rho_mem.astype(float),
             'rho_biomass': 300,
             'phi_mem': mem['mass_frac'].values.astype(float),
             'alpha': size_data['length_um'].values/size_data['width_um'].values,
             'm_peri': ms_data[ms_data['localization'] == 'periplasm']['mass_fg'].values}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)


# %%
# Plot the fit quantites
fig, ax = plt.subplots(2, 2, figsize=(4, 3), sharex=True)
ax = ax.ravel()
ax[0].set_ylim([1, 8])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 6])
ax[3].set_ylim([0, 25])
ax[0].set_ylabel(r'$\alpha$' + '\naspect ratio', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$\nallocation to membrane protein', fontsize=6)
ax[2].set_ylabel(
    r'$\rho_{mem}$ [fg / µm$^2$]' + '\nmembrane protein areal density', fontsize=6)
ax[3].set_ylabel(
    '$m_{peri}$ [fg / cell]\nperiplasmic protein mass', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
for g, d in ms_data[ms_data['localization'].isin(['periplasm', 'membrane'])].groupby(['dataset_name']):
    peri = d[d['localization'] == 'periplasm']
    mem = d[d['localization'] == 'membrane']
    _style = {'markeredgecolor': cor['primary_black'],
              'color': mapper[g]['c'],
              'marker': mapper[g]['m'],
              'linestyle': 'none',
              'label': g,
              'markeredgewidth': 0.25,
              'alpha': 0.5}

    ax[1].plot(mem['growth_rate_hr'], mem['mass_frac'], **_style, zorder=1000)
    ax[2].plot(mem['growth_rate_hr'], mem['mass_fg'] /
               (2 * mem['surface_area']), **_style, zorder=1000)
    ax[3].plot(peri['growth_rate_hr'], peri['mass_fg'], **_style, zorder=1000)

for g, d in size_data.groupby(['source']):
    _style = {'markeredgecolor': cor['primary_black'],
              'color': mapper[g]['c'],
              'marker': mapper[g]['m'],
              'linestyle': 'none',
              'label': g,
              'markeredgewidth': 0.25,
              'alpha': 0.75}
    ax[0].plot(d['growth_rate_hr'], d['length_um'] /
               d['width_um'], **_style, zorder=1000)

# Plot the ppcs
fit_ppc_df = samples.posterior[['alpha_rep', 'm_peri_rep',
                                'rho_mem_rep', 'phi_mem_rep']].to_dataframe(
).reset_index()[['alpha_rep', 'm_peri_rep',
                 'rho_mem_rep', 'phi_mem_rep']
                ].melt()
fit_ppc_df.dropna(inplace=True)
percs = [(2.5, 97.5), (12.5, 87.5), (37.5,  62.5), (50.0, 50.0)]
colors = [cor['pale_black'], cor['light_black'],
          cor['primary_black'], cor['dark_black']]
axes = {'alpha_rep': ax[0],  'phi_mem_rep': ax[1],
        'rho_mem_rep': ax[2], 'm_peri_rep': ax[3]}

for g, d in fit_ppc_df.groupby('variable'):
    for i, p in enumerate(percs):
        lower, upper = np.percentile(d['value'], p)
        axes[g].fill_between(np.linspace(0, 2.5), lower,
                             upper, color=colors[i], alpha=0.25)
plt.tight_layout()
plt.savefig('../../figures/fit_constants.pdf', bbox_inches='tight')

# %%
# Predicted parameters
fig, ax = plt.subplots(2, 2, figsize=(4, 3), sharex=True)
ax = ax.ravel()
ax[0].set_ylim([0.3, 1.5])
ax[1].set_ylim([1, 6])
ax[2].set_ylim([-0.5, 6])
ax[3].set_ylim([0, 150])
ax[0].set_ylabel(f'w [µm]\naverage cell width', fontsize=6)
ax[1].set_ylabel(f'$\ell$ [µm]\naverage cell length', fontsize=6)
ax[2].set_ylabel(f'V [µm$^3$]\naverage cell volume', fontsize=6)
ax[3].set_ylabel(r'$\rho_{peri}$ [fg / µm$^{-3}$]' +
                 '\nperiplasmic protein density', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
for g, d in size_data.groupby('source'):
    _style = {'markeredgecolor': cor['primary_black'],
              'color': mapper[g]['c'],
              'marker': mapper[g]['m'],
              'linestyle': 'none',
              'label': g,
              'markeredgewidth': 0.25,
              'alpha': 0.75}
    ax[0].plot(d['growth_rate_hr'], d['width_um'], **_style)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], **_style)
    ax[2].plot(d['growth_rate_hr'], d['volume_um3'], **_style)

for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby('dataset_name'):
    _style = {'markeredgecolor': cor['primary_black'],
              'color': mapper[g]['c'],
              'marker': mapper[g]['m'],
              'linestyle': 'none',
              'label': g,
              'markeredgewidth': 0.25,
              'alpha': 0.5}
    ax[3].plot(d['growth_rate_hr'], d['mass_fg'].values /
               (d['surface_area'] * 0.0246), **_style)

pred_props = ['width_rep', 'ell_rep', 'vol_rep', 'rho_peri']
percs = [(2.5, 97.5), (12.5, 87.5), (37.5,  62.5), (50.0, 50.0)]
colors = [cor['pale_blue'], cor['light_blue'],
          cor['primary_blue'], cor['dark_blue']]
axes = {'width_rep': ax[0],  'ell_rep': ax[1],
        'vol_rep': ax[2], 'rho_peri': ax[3]}

for i, var in enumerate(pred_props):
    lower = []
    upper = []
    df = samples.posterior[var].to_dataframe().reset_index()
    for j, p in enumerate(percs):
        lower, upper = [], []
        for g, d in df.groupby(f'{var}_dim_0'):
            d.dropna(inplace=True)
            _perc = np.percentile(d[var], p)
            lower.append(_perc[0])
            upper.append(_perc[1])
        axes[var].fill_between(lam_range, lower, upper,
                               color=colors[j], alpha=0.25)

plt.tight_layout()
plt.savefig('../../figures/allocation_pred_ppcs.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 2, figsize=(4, 2))
ax[1].plot(si_clim_data['phi_Rb'], si_clim_data['width_um'], mapper['Si et al., 2017']['m'], label='__nolegend__',
           markeredgecolor=cor['primary_black'], color=mapper['Si et al., 2017']['c'], markeredgewidth=0.25, alpha=0.5)
ax[0].plot(si_clim_data['growth_rate_hr'], si_clim_data['phi_Rb'], mapper['Si et al., 2017']['m'], label='__nolegend__',
           markeredgecolor=cor['primary_black'], color=mapper['Si et al., 2017']['c'], markeredgewidth=0.25, alpha=0.5)

for g, d in ribo_data.groupby('dataset_name'):
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], markerfacecolor=mapper[g]['c'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.25, alpha=0.5, label='__nolegend__')

    ax[1].plot(d['mass_frac'], d['width'], mapper[g]['m'], markerfacecolor=mapper[g]['c'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.25, alpha=0.5, label='__nolegend__')
for g, d in size_data.groupby(['source']):
    if g == 'Si et al., 2017':
        continue
    ax[0].plot(d['growth_rate_hr'], d['phiRb'], mapper[g]['m'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.25, color=mapper[g]['c'], alpha=0.25, label='__nolegend__')

    ax[1].plot(d['phiRb'], d['width_um'], mapper[g]['m'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.25, color=mapper[g]['c'], alpha=0.25, label='__nolegend__')
for g, d in phiRb_data.groupby(['source']):
    if g == 'Si et al. 2017':
        continue
    try:
        ax[0].plot(d['growth_rate_hr'], d['mass_fraction'], markercolors[g]['m'], markeredgecolor=cor['primary_black'],
                   color=markercolors[g]['c'], alpha=0.25, markeredgewidth=0.25, label='__nolegend__')

        ax[1].plot(d['mass_fraction'], d['width'], markercolors[g]['m'], markeredgecolor=cor['primary_black'],
                   color=markercolors[g]['c'], alpha=0.25, markeredgewidth=0.25, label='__nolegend__')
    except:
        continue

ax[0].plot(lam_range, phiRb_range, '-', lw=1,
           color=cor['primary_red'], label='optimal allocation')
df = samples.posterior['width_rep'].to_dataframe().reset_index()
for j, p in enumerate(percs):
    lower, upper = [], []
    for g, d in df.groupby('width_rep_dim_0'):
        d.dropna(inplace=True)
        _perc = np.percentile(d['width_rep'], p)
        lower.append(_perc[0])
        upper.append(_perc[1])
    ax[1].fill_between(phiRb_range, lower, upper,
                       color=colors[j], alpha=0.5)

ax[0].legend(fontsize=6)
ax[1].set_ylim([0.45, 1.2])
ax[1].set_xlim([0, 0.3])
ax[1].set_xlabel('ribosomal mass fraction\n$\phi_{Rb}$', fontsize=6)
ax[1].set_ylabel('w [µm]\naverage cell width', fontsize=6)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('ribosomal mass fraction\n$\phi_{Rb}$', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/allocation_width_prediction.pdf',
            bbox_inches='tight')
# %%
