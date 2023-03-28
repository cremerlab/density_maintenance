# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.analytical
import seaborn as sns
mapper = size.viz.lit_mapper()
err_widths = {'95%': 0.25, '75%': 1, '25%': 2}
cor, pal = size.viz.matplotlib_style()
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
# mass_fracs = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
growth_params = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
params = pd.concat([params, growth_params], sort=False)
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
lines = pd.read_csv('../../data/mcmc/growth_rate_linear_relations.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')
pred_phi = pd.read_csv('../../data/mcmc/predicted_phi_scaling_lppwt.csv')
singular_posts = pd.read_csv('../../data/mcmc/singular_posterior_kde.csv')
singular_percs = pd.read_csv(
    '../../data/mcmc/singular_parameter_percentiles.csv')
mass_spec_medians = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')
mass_spec_medians = mass_spec_medians[mass_spec_medians['interval'] == 'median']

params = params[(params['strain'] == 'wildtype') &
                (params['overexpression'] == 'none') &
                (params['inducer_conc'] == 0)]
posts = posts[(posts['strain'] == 'wildtype') &
              (posts['overexpression'] == 'none') &
              (posts['inducer_conc_ng_mL'] == 0)]
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'})
shapes = shapes[(shapes['strain'] == 'wildtype') &
                (shapes['overexpression'] == 'none') &
                (shapes['inducer_conc'] == 0)]
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']

# %%
# Define the parameters
k = np.array([1, 2, 3])
delta = 0.025
slope = k * delta
phi_range = np.linspace(0.01, 0.10, 100)
sav_range = np.linspace(5.7001, 7.2, 100)

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.75))

ls = ['-', '--', ':', '-.']
for i, s in enumerate(slope):
    theo = s * (sav_range - 5.7)
    ax.plot(sav_range, theo + 0.02,  ls=ls[i],
            lw=1, color=cor['primary_blue'], zorder=1000)

# Plot the mass spec data
for g, d in mass_spec_medians.groupby(['dataset_name']):
    ax.plot(d[d['quantity'] == 'mass_spec_sav']['lower'],
            d[d['quantity'] == 'mass_spec_phi_M']['lower'],
            linestyle='none', marker=mapper[g]['m'], markeredgecolor='k',
            alpha=0.15, markeredgewidth=0.5, color=mapper[g]['c'], ms=4, zorder=10)


med_params = params[(params['quantity'].isin(['surface_area_vol_mu', 'phi_M'])) &
                    (params['interval'] == 'median')]
for g, d in params[params['quantity'].isin(['surface_area_vol_mu', 'phi_M'])].groupby(['carbon_source',
                                                                                       'interval'], sort=False):
    if g[0] == 'LB':
        continue

    if g[1] == 'median':
        ax.plot(d[d['quantity'] == 'surface_area_vol_mu']['lower'], d[d['quantity'] == 'phi_M']['upper'],
                marker='o', markeredgewidth=0.75, markeredgecolor=cor['primary_blue'],
                markerfacecolor='white', ms=3, zorder=1000)
    else:
        sav = d[d['quantity'] == 'surface_area_vol_mu']
        phiM = d[d['quantity'] == 'phi_M']
        med_sav = med_params[(med_params['quantity'] == 'surface_area_vol_mu') &
                             (med_params['carbon_source'] == g[0])]
        med_phi_M = med_params[(med_params['quantity'] == 'phi_M') &
                               (med_params['carbon_source'] == g[0])]
        ax.vlines(med_sav['lower'], phiM['lower'], phiM['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)
        ax.hlines(med_phi_M['lower'], sav['lower'], sav['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)

ax.set_xlim([5.5, 7.2])
ax.set_ylim([0.01, 0.12])
ax.set_xlabel('surface area to volume [µm$^{-1}$]', fontsize=6)
ax.set_ylabel(
    'periplasmic biomass fraction', fontsize=6)
plt.savefig('../../figures/Fig2_sav_scaling.pdf', bbox_inches='tight')


# %%

cors = {'width_slope': 'black', 'length_slope': 'purple', 'alpha_slope': 'blue'}
inters = {'95%': 'pale', '75%': 'light', '25%': 'primary', 'median': ''}
meds = singulars[singulars['interval'] == 'median']
fig, ax = plt.subplots(1, 2, figsize=(3, 1))

for g, d in lit_size_data.groupby(['source']):
    d['aspect_ratio'] = d['length_um'] / d['width_um']
    ax[0].plot(d['growth_rate_hr'], d['aspect_ratio'],
               marker=mapper[g]['m'], ms=3, color=mapper[g]['c'],
               markeredgecolor=cor['primary_black'], alpha=0.3, linestyle='none')

for g, d in params[(params['quantity'].isin(['growth_mu', 'aspect_ratio_mu'])) &
                   (params['strain'] == 'wildtype')].groupby(['carbon_source']):
    med_lam = d[(d['quantity'] == 'growth_mu')
                & (d['interval'] == 'median')]['lower']
    med_alpha = d[(d['quantity'] == 'aspect_ratio_mu')
                  & (d['interval'] == 'median')]['lower']
    ax[0].plot(med_lam, med_alpha, 'o', ms=3, markerfacecolor='white', markeredgecolor=cor['blue'],
               markeredgewidth=0.5, zorder=1000)
    for _g, _d in d[d['interval'] != 'median'].groupby(['interval'], sort=False):
        lam = _d[_d['quantity'] == 'growth_mu']
        alpha = _d[_d['quantity'] == 'aspect_ratio_mu']
        ax[0].hlines(med_alpha, lam['lower'], lam['upper'], lw=err_widths[_g],
                     color=cor['blue'])
        ax[0].vlines(med_lam, alpha['lower'], alpha['upper'], lw=err_widths[_g],
                     color=cor['blue'])

inter_colors = {'95%': cor['pale_blue'],
                '75%': cor['light_blue'], '25%': cor['primary_blue']}
for g, d in lines[(lines['quantity'] == 'alpha') & (lines['interval'] != 'median')
                  ].groupby(['interval'], sort=False):
    ax[0].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                       color=inter_colors[g])

locs = {'alpha': 0, 'width': 1, 'length': 2}
for i, (g, d) in enumerate(singulars[singulars['quantity'].isin(['width_slope', 'length_slope',
                                                                 'alpha_slope'])].groupby(['quantity'])):

    param = g.split('_')[0]
    med_min = singulars[(singulars['quantity'] == f'{param}_min') &
                        (singulars['interval'] == 'median')]['lower'].values[0]
    med = d[d['interval'] == 'median']
    points = med['lower'] / med_min
    ax[1].plot(locs[param], points, 'o', ms=4, markeredgewidth=0.5,
               markeredgecolor=cor['blue'], markerfacecolor='white', zorder=1000)
    for _g, _d in d[d['interval'] != 'median'].groupby(['interval']):
        ax[1].vlines(locs[param], _d['lower']/med_min, _d['upper']/med_min, color=cor['blue'],
                     lw=err_widths[_g])

ax[0].hlines(3.3, 0, 2.5, linestyle='--', color=cor['primary_blue'], lw=1)
ax[0].set_ylim([1, 8])
ax[0].set_ylabel('aspect ratio', fontsize=6)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_xlim([0, 0.90])
ax[1].set_xlim([-1, 3])
ax[1].set_ylabel('relative slope\n[fold-increase per hr]', fontsize=6)
ax[1].set_xticks([0, 1, 2])
plt.savefig('../../figures/Fig2_aspect_ratio_trend.pdf', bbox_inches='tight')

# %%
# Using aspect ratio to compute length and volume
fig, ax = plt.subplots(1, 2, figsize=(3, 1))
for g, d in lit_size_data.groupby(['source']):
    ax[0].plot(d['growth_rate_hr'], d['length_um'], linestyle='none', marker=mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], ms=3, alpha=0.3,
               markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)
    ax[1].plot(d['growth_rate_hr'], d['volume_um3'], linestyle='none', marker=mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], ms=3, alpha=0.3,
               markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5)

singular_meds = singulars[singulars['interval'] == 'median']
lam_range = np.linspace(0, 2.5, 100)
alpha = 3.3
width_trend = singular_meds[singular_meds['quantity'] == 'width_min']['lower'].values[0] + \
    singular_meds[singular_meds['quantity']
                  == 'width_slope']['lower'].values[0] * lam_range
length_trend = alpha * width_trend
vol_trend = (np.pi / 12) * width_trend**2 * (3 * length_trend - width_trend)
ax[0].plot(lam_range, length_trend, '--', color=cor['primary_blue'], lw=1)
ax[1].plot(lam_range, vol_trend, '--', color=cor['primary_blue'], lw=1)

for g, d in params[(params['quantity'].isin(['growth_mu', 'length_mu', 'volume_mu'])) &
                   (params['strain'] == 'wildtype')].groupby(['carbon_source']):
    med_lam = d[(d['quantity'] == 'growth_mu')
                & (d['interval'] == 'median')]['lower']
    med_ell = d[(d['quantity'] == 'length_mu')
                & (d['interval'] == 'median')]['lower']
    med_vol = d[(d['quantity'] == 'volume_mu')
                & (d['interval'] == 'median')]['lower']
    ax[0].plot(med_lam, med_ell, 'o', ms=3, markerfacecolor='white', markeredgecolor=cor['blue'],
               markeredgewidth=0.5, zorder=1000)
    ax[1].plot(med_lam, med_vol, 'o', ms=3, markerfacecolor='white', markeredgecolor=cor['blue'],
               markeredgewidth=0.5, zorder=1000)
    for _g, _d in d[d['interval'] != 'median'].groupby(['interval'], sort=False):
        lam = _d[_d['quantity'] == 'growth_mu']
        ell = _d[_d['quantity'] == 'length_mu']
        vol = _d[_d['quantity'] == 'volume_mu']
        ax[0].hlines(med_ell, lam['lower'], lam['upper'], lw=err_widths[_g],
                     color=cor['blue'])
        ax[1].hlines(med_vol, lam['lower'], lam['upper'], lw=err_widths[_g],
                     color=cor['blue'])

        ax[0].vlines(med_lam, ell['lower'], ell['upper'], lw=err_widths[_g],
                     color=cor['blue'])
        ax[1].vlines(med_lam, vol['lower'], vol['upper'], lw=err_widths[_g],
                     color=cor['blue'])

ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('average length [µm]', fontsize=6)
ax[1].set_ylabel('average volume [µm$^{3}$]', fontsize=6)
ax[1].set_ylim([0, 5])
ax[0].set_ylim([0, 7])
ax[0].set_yticks([0, 2, 4, 6])
ax[1].set_yticks([1, 3, 5])
plt.savefig('../../figures/Fig2_aspect_ratio_length_vol.pdf',
            bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 1, figsize=(2.25, 1.75))
k = np.array([1, 2, 3])
delta = 0.025
alpha = 3.3
slope = 12 * alpha * k * delta / (3 * alpha - 1)
width_range = np.linspace(0.6, 1.0, 100)
phi_range = np.linspace(0.01, 0.10, 100)


ls = ['-', '--', ':', '-.']
for i, s in enumerate(slope):
    theo = 0.1 - s * (0.6**-1 - width_range**-1)
    ax.plot(width_range, theo,  ls=ls[i],
            lw=1, color=cor['primary_blue'], zorder=1000)

# Plot the mass spec data
for g, d in mass_spec_medians.groupby(['dataset_name']):
    ax.plot(d[d['quantity'] == 'mass_spec_widths']['lower'],
            d[d['quantity'] == 'mass_spec_phi_M']['lower'],
            linestyle='none', marker=mapper[g]['m'], markeredgecolor='k',
            alpha=0.15, markeredgewidth=0.5, color=mapper[g]['c'], ms=4, zorder=10)


med_params = params[(params['quantity'].isin(['width_mu', 'phi_M'])) &
                    (params['interval'] == 'median')]
for g, d in params[params['quantity'].isin(['width_mu', 'phi_M'])].groupby(['carbon_source',
                                                                            'interval'], sort=False):
    if g[0] == 'LB':
        continue

    if g[1] == 'median':
        ax.plot(d[d['quantity'] == 'width_mu']['lower'], d[d['quantity'] == 'phi_M']['upper'],
                marker='o', markeredgewidth=0.75, markeredgecolor=cor['primary_blue'],
                markerfacecolor='white', ms=3, zorder=1000)
    else:
        width = d[d['quantity'] == 'width_mu']
        phiM = d[d['quantity'] == 'phi_M']
        med_width = med_params[(med_params['quantity'] == 'width_mu') &
                               (med_params['carbon_source'] == g[0])]
        med_phi_M = med_params[(med_params['quantity'] == 'phi_M') &
                               (med_params['carbon_source'] == g[0])]
        ax.vlines(med_width['lower'], phiM['lower'], phiM['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)
        ax.hlines(med_phi_M['lower'], width['lower'], width['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)

ax.set_xlim([0.6, 1])
ax.set_ylim([0.005, 0.10])
ax.set_xlabel('average width [µm$^{-1}$]', fontsize=6)
ax.set_ylabel('periplasmic biomass fraction', fontsize=6)
plt.savefig('../../figures/Fig2_width_scaling.pdf', bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 2))
ax[0].set_ylim([0, 250])
ax[1].set_ylim([0, 50])
meds = params[params['interval'] == 'median']
for g, d in mass_spec_medians.groupby(['dataset_name']):
    rho_peri = d[d['quantity'] == 'mass_spec_rho_peri']
    peri_prot = d[d['quantity'] == 'mass_spec_peri_prot_per_cell']
    ax[0].plot(rho_peri['growth_rate_hr'], rho_peri['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor='k',
               markeredgewidth=0.5, alpha=0.3)
    ax[1].plot(peri_prot['growth_rate_hr'], peri_prot['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5, alpha=0.3)

for g, d in params[params['interval'] != 'median'].groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    _lam = meds[(meds['carbon_source'] == g[0]) &
                (meds['quantity'] == 'growth_mu')]
    _rho_peri = meds[(meds['carbon_source'] == g[0]) &
                     (meds['quantity'] == 'rho_peri')]
    _peri_prot = meds[(meds['carbon_source'] == g[0]) &
                      (meds['quantity'] == 'peri_prot_per_cell')]

    ax[0].plot(_lam['lower'], _rho_peri['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[0].hlines(_rho_peri['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[0].vlines(_lam['lower'], d[d['quantity'] == 'rho_peri']['lower'] * 1E9, d[d['quantity'] == 'rho_peri']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])

    ax[1].plot(_lam['lower'], _peri_prot['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[1].hlines(_peri_prot['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[1].vlines(_lam['lower'], d[d['quantity'] == 'peri_prot_per_cell']['lower'] * 1E9, d[d['quantity'] == 'peri_prot_per_cell']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])


ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('periplasmic protein density [fg / µm$^3$]', fontsize=6)
ax[1].set_ylabel('periplasmic protein mass [fg]', fontsize=6)
plt.savefig('/Users/gchure/Desktop/plots.pdf')
# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 2))
meds = params[params['interval'] == 'median']
ax[0].set_ylim([0, 300])
ax[1].set_ylim([0, 300])
# ax[1].set_xlim([0, 1])
for g, d in mass_spec_medians.groupby(['dataset_name']):
    rho_peri = d[d['quantity'] == 'mass_spec_rho_peri']
    widths = d[d['quantity'] == 'mass_spec_widths']
    ax[0].plot(rho_peri['growth_rate_hr'], rho_peri['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor='k',
               markeredgewidth=0.5, alpha=0.3)
    ax[1].plot(widths['lower'], rho_peri['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5, alpha=0.3)
width_range = np.linspace(0.6, 1)
slope1 = 30 / (np.pi * 3.3 * 0.024)
slope2 = 20 / (np.pi * 3.3 * 0.024)
slope3 = 10 / (np.pi * 3.3 * 0.024)
theo1 = slope1 / width_range**2
theo2 = slope2 / width_range**2
theo3 = slope3 / width_range**2
ax[1].plot(width_range, theo1, 'k-')
ax[1].plot(width_range, theo2, 'k--')
ax[1].plot(width_range, theo3, 'k:')
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('periplasmic protein density [fg / µm$^3$]', fontsize=6)
ax[1].set_xlabel('cell width [µm]', fontsize=6)
ax[1].set_ylabel('periplasmic protein density [fg / µm$^3$]', fontsize=6)
plt.savefig('/Users/gchure/Desktop/scaling_argument.pdf', bbox_inches='tight')
# %%
for g, d in params[params['interval'] != 'median'].groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    _lam = meds[(meds['carbon_source'] == g[0]) &
                (meds['quantity'] == 'growth_mu')]
    _rho_peri = meds[(meds['carbon_source'] == g[0]) &
                     (meds['quantity'] == 'rho_peri')]
    _peri_prot = meds[(meds['carbon_source'] == g[0]) &
                      (meds['quantity'] == 'peri_prot_per_cell')]

    ax[0].plot(_lam['lower'], _rho_peri['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[0].hlines(_rho_peri['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[0].vlines(_lam['lower'], d[d['quantity'] == 'rho_peri']['lower'] * 1E9, d[d['quantity'] == 'rho_peri']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])

    ax[1].plot(_lam['lower'], _peri_prot['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[1].hlines(_peri_prot['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[1].vlines(_lam['lower'], d[d['quantity'] == 'peri_prot_per_cell']['lower'] * 1E9, d[d['quantity'] == 'peri_prot_per_cell']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])


# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 2))
ax[0].set_ylim([0, 250])
ax[1].set_ylim([0, 50])
meds = params[params['interval'] == 'median']
for g, d in mass_spec_medians.groupby(['dataset_name']):
    rho_peri = d[d['quantity'] == 'mass_spec_rho_peri']
    peri_prot = d[d['quantity'] == 'mass_spec_peri_prot_per_cell']
    ax[0].plot(rho_peri['growth_rate_hr'], rho_peri['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor='k',
               markeredgewidth=0.5, alpha=0.3)
    ax[1].plot(peri_prot['growth_rate_hr'], peri_prot['lower'] * 1E9, mapper[g]['m'],
               ms=4, markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5, alpha=0.3)

for g, d in params[params['interval'] != 'median'].groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    _lam = meds[(meds['carbon_source'] == g[0]) &
                (meds['quantity'] == 'growth_mu')]
    _rho_peri = meds[(meds['carbon_source'] == g[0]) &
                     (meds['quantity'] == 'rho_peri')]
    _peri_prot = meds[(meds['carbon_source'] == g[0]) &
                      (meds['quantity'] == 'peri_prot_per_cell')]

    ax[0].plot(_lam['lower'], _rho_peri['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[0].hlines(_rho_peri['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[0].vlines(_lam['lower'], d[d['quantity'] == 'rho_peri']['lower'] * 1E9, d[d['quantity'] == 'rho_peri']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])

    ax[1].plot(_lam['lower'], _peri_prot['lower'] * 1E9, 'o',
               markeredgecolor=cor['primary_blue'], markerfacecolor='white',
               ms=3, markeredgewidth=1, zorder=1000)
    ax[1].hlines(_peri_prot['lower'] * 1E9, d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'growth_mu']['upper'],
                 color=cor['primary_blue'], lw=err_widths[g[1]])
    ax[1].vlines(_lam['lower'], d[d['quantity'] == 'peri_prot_per_cell']['lower'] * 1E9, d[d['quantity'] == 'peri_prot_per_cell']['upper'] * 1E9,
                 color=cor['primary_blue'], lw=err_widths[g[1]])


ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('periplasmic protein density [fg / µm$^3$]', fontsize=6)
ax[1].set_ylabel('periplasmic protein mass [fg]', fontsize=6)
plt.savefig('/Users/gchure/Desktop/plots.pdf')
