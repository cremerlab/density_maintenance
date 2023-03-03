# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import Blues_7
import size.viz
import size.analytical
cor, pal = size.viz.matplotlib_style()

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')

# Restrict to wildtype
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

# Generate a KDE plot of the posteriors for aspect ratio and density ratio
_colors = Blues_7.mpl_colors[2:]
_colors = {c: _colors[i] for i, c in enumerate(
    ['acetate', 'sorbitol', 'glycerol', 'glucose', 'glucoseCAA'])}

# Set up the canvas and label
fig, ax = plt.subplots(1, 2, figsize=(2, 0.75))
for a in ax:
    a.set_yticks([])
    a.xaxis.set_label_position('top')
    a.xaxis.tick_top()
ax[0].set_xlim([2.5, 4])
ax[0].set_xticks([2.5, 3, 3.5])
ax[0].set_xticklabels(['', '', ''])
ax[1].set_xlim([0, 0.4])
ax[1].set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax[1].set_xticklabels(['', '', '', '', ''])


axes = {'aspect_ratio_mu': 0, 'rho_ratio': 1}
for i, (g, d) in enumerate(posts[posts['parameter'].isin(['aspect_ratio_mu', 'rho_ratio'])].groupby(['carbon_source'])):
    if g == 'LB':
        continue
    for _g, _d in d.groupby(['parameter']):
        ax[axes[_g]].fill_between(
            _d['value'], _d['kde'], color=_colors[g], alpha=0.25, zorder=i+1)
        ax[axes[_g]].plot(_d['value'], _d['kde'], '-',
                          color=_colors[g], lw=0.75, zorder=i+1)
plt.tight_layout()
plt.savefig('../../figures/Fig2_constant_parameters_kde.pdf')
# %%
fig, ax = plt.subplots(1, 2, figsize=(2, 1), sharey=True)
ycoords = {'glucoseCAA': 5, 'glucose': 4,
           'glycerol': 3, 'sorbitol': 2, 'acetate': 1}
axes = {'aspect_ratio_mu': 0, 'rho_ratio': 1}
err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}

for g, d in errs[errs['quantity'].isin(list(axes.keys()))
                 ].groupby(['carbon_source', 'quantity', 'interval']):
    if g[0] == 'LB':
        continue
    ax[axes[g[1]]].hlines(ycoords[g[0]], d['lower'], d['upper'], lw=err_widths[g[-1]],
                          color=_colors[g[0]])


for g, d in medians[medians['quantity'].isin(list(axes.keys()))
                    ].groupby(['carbon_source', 'quantity']):
    if g[0] == 'LB':
        continue
    ax[axes[g[1]]].plot(d['lower'], ycoords[g[0]], 'o', markeredgewidth=0.5, ms=3,
                        markeredgecolor=_colors[g[0]], markerfacecolor='white')

# Add labels
ax[0].set_yticks(list(ycoords.values()))
ax[0].set_ylim([0.5, 5.5])
labels = ['glucose\n+ CAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']

# Update limits
# ax[0].set_yticklabels(labels)
ax[0].set_xlim([2.5, 4])
ax[0].set_xticks([2.5, 3, 3.5])
ax[1].set_xlim([0, 0.4])
ax[1].set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax[0].set_xticklabels(['', '', ''])
ax[1].set_xticklabels(['', '', '', '', '', ])
ax[0].set_yticks([])
for a in ax:
    a.xaxis.tick_top()
# # Add labels
plt.tight_layout()
plt.savefig('../../figures/Fig2_constant_parameters_percentiles.pdf',
            bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.5), sharey=True)

# Plot the percentiles and the medians
for g, d in errs[(errs['quantity'].isin(['phi_M', 'width_mu']))
                 ].groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    phi = medians[(medians['carbon_source'] == g[0]) &
                  (medians['quantity'] == 'phi_M')]['lower']
    w = medians[(medians['carbon_source'] == g[0]) & (
        medians['quantity'] == 'width_mu')]['lower']
    phi_d = d[d['quantity'] == 'phi_M']
    w_d = d[d['quantity'] == 'width_mu']
    ax.hlines(w, phi_d['lower'], phi_d['upper'], lw=err_widths[g[1]],
              color=cor['primary_blue'])
    ax.vlines(phi, w_d['lower'], w_d['upper'], lw=err_widths[g[1]],
              color=cor['primary_blue'])

for g, d in medians[medians['quantity'].isin(['phi_M',
                                              'width_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    phi = d[d['quantity'] == 'phi_M']['lower']
    w = d[d['quantity'] == 'width_mu']['lower']
    ax.plot(phi, w, 'o', markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
            markerfacecolor='white', ms=3)
# phi_range = np.linspace(0, 0.06, 100)
# theo = 12 * 3 * 0.024 * 0.17 / (8 * phi_range)
# ax.plot(phi_range, theo)
# ax.set_ylim([0.45, 1.2])
band_cor = {'95%': cor['pale_black'],
            '75%': cor['light_black'], '25%': cor['primary_black']}
for i, (g, d) in enumerate(pred_err.groupby(['interval'], sort=False)):
    ax.fill_between(d['phi_M'], d['lower'], d['upper'],
                    color=band_cor[g], alpha=0.75)
ax.set_ylim([0.45, 1])
ax.set_xlim([0, 0.06])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_ylabel('w\naverage width [µm]', fontsize=6)
plt.savefig('../../figures/Fig2_width_prediction_wildtype.pdf')
# %%
phi_range = np.linspace(0, 0.06, 10)
k = 0.17
delta = 0.024
theo = (1 + (k / (phi_range**-1 - 1)))**-1
plt.plot(phi_range, theo)
# plt.ylim([0, 10])

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))

phi_range = np.linspace(0, 0.06, 100)
x = [0.1, 0.15, 0.2]
delta = 0.024
for i, _x in enumerate(x):
    if i == 0:
        ls = '--'
    elif i == 1:
        ls = '-.'
    else:
        ls = ':'
    theo = (delta * (_x * (1/phi_range - 1) + 1))**-1
    ax.plot(phi_range, theo, ls=ls, lw=1, color=cor['light_black'])

for g, d in errs.groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    phi_loc = medians[(medians['quantity'] == 'phi_M') & (
        medians['carbon_source'] == g[0])]['lower']
    sv_loc = medians[(medians['quantity'] == 'surface_area_vol_mu') &
                     (medians['carbon_source'] == g[0])]['lower']
    phi = d[d['quantity'] == 'phi_M']
    sv = d[d['quantity'] == 'surface_area_vol_mu']
    ax.hlines(sv_loc, phi['lower'], phi['upper'], color=cor['primary_blue'],
              lw=err_widths[g[1]])
    ax.vlines(phi_loc, sv['lower'], sv['upper'], color=cor['primary_blue'],
              lw=err_widths[g[1]])

for g, d in medians.groupby(['carbon_source']):
    if g == 'LB':
        continue
    ax.plot(d[d['quantity'] == 'phi_M']['lower'],
            d[d['quantity'] == 'surface_area_vol_mu']['lower'], 'o',
            ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
            markerfacecolor='w')


ax.set_ylim([4, 8])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_{M}$', fontsize=6)
ax.set_ylabel('$S/V$\nsurface area to volume [µm$^{-1}$]', fontsize=6)
plt.savefig('../../figures/Fig2_SV_scaling.pdf', bbox_inches='tight')
