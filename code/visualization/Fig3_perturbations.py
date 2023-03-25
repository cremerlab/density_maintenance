# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import Blues_7
import size.viz
import size.analytical
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}
# Load datasets
growth_params = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
params = pd.concat([params, growth_params], sort=False)
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
mass_spec_medians = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')
mass_spec_medians = mass_spec_medians[mass_spec_medians['interval'] == 'median']


# Restrict to wildtype
params = params[~((params['overexpression'] != 'none')
                  * (params['overexpression'] == 0))]
posts = posts[~((posts['overexpression'] != 'none')
                * (posts['overexpression'] == 0))]
params.loc[params['quantity'] == 'prot_per_cell', 'lower'] *= (1E-6/1E-15)
params.loc[params['quantity'] == 'prot_per_cell', 'upper'] *= (1E-6/1E-15)
posts.loc[posts['parameter'] == 'prot_per_cell', 'value'] *= (1E-6/1E-15)
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'}, inplace=True)
shapes = shapes[~((shapes['overexpression'] != 'none')
                  * (shapes['overexpression'] == 0))]
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']

# %%
fig, ax = plt.subplots(1, 3, figsize=(3, 2))
for a in ax.ravel():
    a.set_yticks([])
ax[0].set_xlim([0.01, 0.06])
ax[1].set_xlim([0.5, 0.9])
ax[2].set_xlim([2, 4])
ax[0].set_xticks([0.02, 0.04, 0.06])
ax[1].set_xticks([0.6, 0.7, 0.8])
ax[0].set_xlabel('periplasmic\nbiomass fraction', fontsize=6)
ax[1].set_xlabel('average width\n[µm]', fontsize=6)
ax[2].set_xlabel('aspect ratio', fontsize=6)
n_rows = 5
for i, carb in enumerate(['acetate']):
    for j, param in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
        for k, strain in enumerate(['wildtype', 'malE-rbsB-fliC-KO']):
            samp = posts[(posts['strain'] == strain) &
                         (posts['carbon_source'] == carb) &
                         (posts['overexpression'] == 'none') &
                         (posts['inducer_conc'] == 0) &
                         (posts['parameter'] == param)
                         ]
            percs = params[(params['strain'] == strain) &
                           (params['carbon_source'] == carb) &
                           (params['overexpression'] == 'none') &
                           (params['inducer_conc'] == 0) &
                           (params['quantity'] == param)]
            if strain == 'wildtype':
                for ell in range(n_rows):
                    if ell == 0:
                        ls = '-'
                        alpha = 0.5
                    else:
                        ls = '--'
                        alpha = 0

                    ax[j].fill_between(samp['value'], (n_rows - ell) * np.ones(len(samp)), (n_rows - ell) +
                                       samp['kde'] / samp['kde'].max(),
                                       color=cor['pale_blue'],
                                       alpha=alpha)
                    ax[j].plot(samp['value'], (n_rows - ell) + samp['kde'] / samp['kde'].max(),
                               ls=ls, color=cor['primary_blue'], lw=0.5)
                    for g, d in percs.groupby(['interval'], sort=False):
                        if g != 'median':
                            ax[j].hlines((n_rows - ell) + 0.5, d['lower'], d['upper'],
                                         lw=err_widths[g], color=cor['primary_blue'], zorder=999)
                        else:
                            ax[j].plot(d['lower'], (n_rows - ell) + 0.5, 'o',
                                       ms=3, markeredgecolor=cor['primary_blue'],
                                       markerfacecolor='white', markeredgewidth=0.5, zorder=1000)

            else:
                colors = {'none': 'purple', 'malE': 'green'}
                for oe in ['none', 'malE']:
                    if oe == 'none':
                        ax[j].fill_between(samp['value'], (n_rows - 1) * np.ones(len(samp)), (n_rows - 1) +
                                           samp['kde'] / samp['kde'].max(),
                                           color=cor[f'pale_{colors[oe]}'],
                                           alpha=0.5)
                        ax[j].plot(samp['value'], (n_rows - 1) + samp['kde'] / samp['kde'].max(),
                                   color=cor[f'primary_{colors[oe]}'], lw=0.5)
                        for g, d in percs.groupby(['interval'], sort=False):
                            if g != 'median':
                                ax[j].hlines((n_rows - 1) + 0.3, d['lower'], d['upper'],
                                             lw=err_widths[g], color=cor[f'primary_{colors[oe]}'], zorder=999)
                            else:
                                ax[j].plot(d['lower'], (n_rows - 1) + 0.3, 'o',
                                           ms=3, markeredgecolor=cor[f'primary_{colors[oe]}'],
                                           markerfacecolor='white', markeredgewidth=0.5, zorder=1000)

                    else:
                        samp = posts[(posts['strain'] == strain) &
                                     (posts['carbon_source'] == carb) &
                                     (posts['overexpression'] == oe) &
                                     (posts['parameter'] == param)
                                     ]
                        percs = params[(params['strain'] == strain) &
                                       (params['carbon_source'] == carb) &
                                       (params['overexpression'] == oe) &
                                       (params['quantity'] == param)]
                        for ell, (g, d) in enumerate(samp.groupby(['inducer_conc'])):
                            print(g)
                            ax[j].fill_between(d['value'], (n_rows - 2 - ell) * np.ones(len(d)), (n_rows - 2 - ell) +
                                               d['kde'] /
                                               d['kde'].max(),
                                               color=cor[f'pale_{colors[oe]}'],
                                               alpha=0.5)
                            ax[j].plot(d['value'], (n_rows - 2 - ell) + d['kde'] / d['kde'].max(),
                                       color=cor[f'primary_{colors[oe]}'], lw=0.5)
                            for _g, _d in percs[percs['inducer_conc'] == g].groupby(['interval'], sort=False):
                                if _g != 'median':
                                    ax[j].hlines((n_rows - 2 - ell) + 0.3, _d['lower'], _d['upper'],
                                                 lw=err_widths[_g], color=cor[f'primary_{colors[oe]}'], zorder=999)
                                else:
                                    ax[j].plot(_d['lower'], (n_rows - 2 - ell) + 0.3, 'o',
                                               ms=3, markeredgecolor=cor[f'primary_{colors[oe]}'],
                                               markerfacecolor='white', markeredgewidth=0.5, zorder=1000)
# plt.savefig('../../figures/Fig3_example_posterior_migration.pdf')
# %%

fig, ax = plt.subplots(1, 1, figsize=(2.25, 2))
k = np.array([1, 2, 3])
delta = 0.025
alpha = 3.3
slope = 12 * alpha * k * delta / (3 * alpha - 1)
width_range = np.linspace(0.6, 1.0, 100)
phi_range = np.linspace(0.01, 0.10, 100)

colors = {'none': 'purple', 'malE': 'green', 'rbsB': 'gold'}
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
for g, d in params[(params['quantity'].isin(['width_mu', 'phi_M'])) &
                   (params['overexpression'].isin(['none', 'malE']))].groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source',
                                                                               'interval'], sort=False):
    if g[0] == 'wildtype':
        c = cor['primary_blue']
    else:
        c = cor[f'primary_{colors[g[1]]}']
    if g[3] == 'LB':
        continue

    if g[4] == 'median':
        ax.plot(d[d['quantity'] == 'width_mu']['lower'], d[d['quantity'] == 'phi_M']['upper'],
                marker='o', markeredgewidth=0.75, markeredgecolor=c,
                markerfacecolor='white', ms=3, zorder=1000)
    else:
        width = d[d['quantity'] == 'width_mu']
        phiM = d[d['quantity'] == 'phi_M']
        med_width = med_params[(med_params['quantity'] == 'width_mu') &
                               (med_params['strain'] == g[0]) &
                               (med_params['overexpression'] == g[1]) &
                               (med_params['inducer_conc'] == g[2]) &
                               (med_params['carbon_source'] == g[3])]
        med_phi_M = med_params[(med_params['quantity'] == 'phi_M') &
                               (med_params['strain'] == g[0]) &
                               (med_params['overexpression'] == g[1]) &
                               (med_params['inducer_conc'] == g[2]) &
                               (med_params['carbon_source'] == g[3])]
        ax.vlines(med_width['lower'], phiM['lower'], phiM['upper'], lw=err_widths[g[-1]],
                  color=c, zorder=99)
        ax.hlines(med_phi_M['lower'], width['lower'], width['upper'], lw=err_widths[g[-1]],
                  color=c, zorder=99)

ax.set_xlim([0.6, 1])
ax.set_ylim([0.005, 0.10])
ax.set_xlabel('average width [µm$^{-1}$]', fontsize=6)
ax.set_ylabel('periplasmic biomass fraction', fontsize=6)
plt.savefig('../../figures/Fig3_width_scaling_overexpression.pdf',
            bbox_inches='tight')

# %%

fig, ax = plt.subplots(2, 2, figsize=(3, 3))

ax = ax.ravel()
ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 6])
ax[2].set_ylim([1, 8])
# Plot the mass spec data
for g, d in lit_size_data.groupby(['source']):
    ax[0].plot(d['growth_rate_hr'], d['width_um'], mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               alpha=0.35, ms=3)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               alpha=0.35, ms=3)
    ax[2].plot(d['growth_rate_hr'], d['length_um'] / d['width_um'], mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               alpha=0.35, ms=3)
    ax[3].plot(d['growth_rate_hr'], d['surface_to_volume'], mapper[g]['m'],
               markerfacecolor=mapper[g]['c'], markeredgecolor=cor['primary_black'],
               alpha=0.35, ms=3)

for g, d in params[(params['quantity'].isin(['surface_area_vol_mu', 'growth_mu', 'width_mu', 'length_mu', 'aspect_ratio_mu'])) &
                   (params['overexpression'].isin(['none', 'malE']))].groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source',
                                                                               'interval'], sort=False):
    if g[0] == 'wildtype':
        c = cor['primary_blue']
    else:
        c = cor[f'primary_{colors[g[1]]}']
    # if g[3] == 'LB':
        # continue

    if g[4] == 'median':

        ax[0].plot(d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'width_mu']['upper'],
                   marker='o', markeredgewidth=0.75, markeredgecolor=c,
                   markerfacecolor='white', ms=3, zorder=1000)
        ax[1].plot(d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'length_mu']['upper'],
                   marker='o', markeredgewidth=0.75, markeredgecolor=c,
                   markerfacecolor='white', ms=3, zorder=1000)
        ax[2].plot(d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'aspect_ratio_mu']['upper'],
                   marker='o', markeredgewidth=0.75, markeredgecolor=c,
                   markerfacecolor='white', ms=3, zorder=1000)
        ax[3].plot(d[d['quantity'] == 'growth_mu']['lower'], d[d['quantity'] == 'surface_area_vol_mu']['upper'],
                   marker='o', markeredgewidth=0.75, markeredgecolor=c,
                   markerfacecolor='white', ms=3, zorder=1000)

# %%
