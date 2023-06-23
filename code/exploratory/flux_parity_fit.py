# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import size.fluxparity
import scipy.stats
import arviz as az
mapper = size.viz.lit_mapper()
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
biomass_density = pd.read_csv(
    '../../data/literature/collated_drymass_densities.csv')
ms_data.drop_duplicates(inplace=True)
# %%
width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
size_data['mass_fraction'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width_um'] = width_popt[1] + \
    width_popt[0] * phiRb_data['growth_rate_hr']
ms = ms_data[ms_data['localization'] == 'ribosomal sector'].copy()
ms.rename(columns={'dataset_name': 'source', 'width': 'width_um',
          'mass_frac': 'mass_fraction'}, inplace=True)
merged = pd.concat([size_data[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']],
                   phiRb_data[['source', 'width_um',
                               'mass_fraction', 'growth_rate_hr']],
                   ms[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']]])

# %%
mem = ms_data[ms_data['localization'] == 'membrane'].drop_duplicates()
rho_mem = mem['mass_fg'] / (2 * mem['surface_area'])
N_ppc = 500
nu_range = np.linspace(0.1, 20, N_ppc)
phiRb_range = size.fluxparity.phiRb_optimal_allocation(consts['gamma_max'], nu_range,
                                                       consts['Kd_cpc'], consts['phi_O'])
lam_range = size.fluxparity.steady_state_growth_rate(
    consts['gamma_max'], phiRb_range, nu_range, consts['Kd_cpc'], consts['phi_O'])
# lam_range = np.linspace(0, 2.5, 500)
# phiRb_range = phiRb_popt[1] + phiRb_popt[0] * lam_range

model = cmdstanpy.CmdStanModel(
    stan_file="membrane_density_scaling_inference.stan")
# %%
data_dict = {'N_ms': len(mem),
             'N_size': len(size_data),
             'N_biomass': len(biomass_density),
             'N_phiRb':  len(merged),
             'phiRb_lam': merged['growth_rate_hr'].values,
             'N_ppc': N_ppc,
             'phiRb_range': phiRb_range,
             'delta': 0.0246,
             'beta_rp': 0.4558,
             'rho_biomass': biomass_density['drymass_density_fg_fL'].values,
             'alpha': size_data['length_um'].values/size_data['width_um'].values,
             'widths': merged['width_um'].values,
             'phiRb': merged['mass_fraction'].values,
             'size_lam': size_data['growth_rate_hr'].values,
             'linear_alpha': 1,
             'lam_range': lam_range,
             'm_peri': ms_data[ms_data['localization'] == 'periplasm']['mass_fg'].values}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(2, 4, figsize=(8, 4))
ax = ax.ravel()
ax[0].axis('off')
_ax = ax[0]
ax = ax[1:]
ax[0].set_xlabel('allocation towards ribosomes\n$\phi_{Rb}$', fontsize=6)
ax[0].set_ylabel('w [µm]\naverage cell width', fontsize=6)
ax[3].set_ylabel(r'$\ell$ [µm]' + '\naverage cell length', fontsize=6)
ax[4].set_ylabel(r'$V$ [µm$^{3}$]' + '\naverage cell volume', fontsize=6)
ax[5].set_ylabel(r'$\rho_{peri}$ [fg / µm$^3$]' +
                 '\nperiplasmic protein density', fontsize=6)
for a in ax[1:]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

# axis  limits
ax[0].set_ylim([0.4, 1.2])
ax[0].set_xlim([0, 0.3])
ax[3].set_ylim([1, 6])
ax[4].set_xlim([0, 2.65])
ax[5].set_ylim([-0.5, 6])
ax[4].set_xlim([0, 2.65])
ax[5].set_ylim([0, 150])

for g, d in merged.groupby('source'):
    _style = size.viz.style_point(g)
    ax[0].plot(d['mass_fraction'], d['width_um'], **_style)

for g, d in size_data.groupby('source'):
    _style = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['length_um'] / d['width_um'], **_style)
    ax[3].plot(d['growth_rate_hr'], d['length_um'], **_style)
    ax[4].plot(d['growth_rate_hr'], d['volume_um3'], **_style)

for g, d in ms_data.groupby('dataset_name'):
    _style = size.viz.style_point(g)
    im_style = _style.copy()
    om_style = _style.copy()
    im_style['color'] = cor['light_red']
    om_style['color'] = cor['light_purple']
    peri = d[d['localization'] == 'periplasm']
    im = d[d['localization'] == 'inner membrane']
    im['zeta'] = im['mass_fg'] / (im['surface_area'] * im['mass_frac'])
    om = d[d['localization'] == 'outer membrane']
    om['zeta'] = om['mass_fg'] / (om['surface_area'] * om['mass_frac'])
    mem = d[d['localization'] == 'membrane']
    mem['zeta'] = mem['mass_fg'] / (2 * mem['surface_area'] * mem['mass_frac'])
    ax[5].plot(peri['growth_rate_hr'], peri['mass_fg'] /
               (0.0246 * peri['surface_area']), **_style)
    # ax[6].plot(im['growth_rate_hr'], im['zeta'], **im_style)
    ax[6].plot(om['growth_rate_hr'], om['zeta'], **om_style)
    # ax[6].plot(mem['growth_rate_hr'], mem['zeta'], **_style)


percs = [(2.5, 97.5), (12.5, 87.5), (37.5,  62.5), (50.0, 50.0)]
hue = ['pale', 'light', 'primary', 'dark']
color = ['black', 'blue', 'blue', 'blue']
ppcs = ['w_rep', 'ell_rep', 'vol_rep', 'rho_peri_rep']
axes = {'w_rep': ax[0], 'ell_rep': ax[3],
        'vol_rep': ax[4], 'rho_peri_rep': ax[5]}

for i, par in enumerate(ppcs):
    df = samples.posterior[par].to_dataframe().reset_index()
    for j, p in enumerate(percs):
        lower, upper = [], []
        for g, d in df.groupby(f'{par}_dim_0'):
            _p = np.percentile(d[par], p)
            lower.append(_p[0])
            upper.append(_p[1])
        if i == 0:
            x = phiRb_range
        else:
            x = lam_range
        axes[par].fill_between(x, lower, upper, color=cor[f'{hue[j]}_{color[i]}'],
                               alpha=0.5)

plt.tight_layout()

# %%
