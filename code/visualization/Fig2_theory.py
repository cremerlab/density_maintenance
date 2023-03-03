# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import Blues_7
import size.viz
import size.analytical
cor, pal = size.viz.matplotlib_style()

# %%

delta = 0.029
width_range = np.linspace(0.5, 3, 200)
ref_width = 1
ref_length = 3
ref_volume = (np.pi/12) * ref_width**2 * (3 * ref_length - ref_width)
volume = (np.pi/12) * width_range**2 * (3 * ref_length - ref_width)
peri_vol = delta * np.pi * ref_length * width_range
ref_peri_vol = delta * np.pi * ref_length * ref_width
plt.plot(width_range, volume/ref_volume, color='blue')
plt.plot(width_range, peri_vol/ref_peri_vol, color='purple')

# %%
length_range = np.linspace(0.5, 4, 200)
ref_width = 1
ref_length = 3
width_range = np.linspace(0.25 * ref_width, 2 * ref_width, 200)
length_range = np.linspace(ref_width, 2 * ref_length, 200)
ref_SAV = 12 * ref_length / (ref_width * (3 * ref_length - ref_width))
width_SAV = 12 * ref_length / (width_range * (3 * ref_length - width_range))
length_SAV = 12 * length_range / (ref_width * (3 * length_range - ref_width))

fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.plot(width_range/ref_width, width_SAV /
        ref_SAV, color=cor['primary_blue'], lw=1)
ax.plot(length_range/ref_length, length_SAV /
        ref_SAV, color=cor['primary_purple'], lw=1)
ax.set_ylim([0.5, 2])
ax.set_xlim([0.3, 1.75])
ax.set_yticks([0.5, 1.0, 1.5, 2.0])
plt.savefig('../../figures/Fig2_length_width_SAV_sensitivity.pdf')

# %%
delta = 0.029
# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')

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
                (shapes['inducer_conc'] == 0) &
                (shapes['parameter'] == 'aspect_ratio_mu')]

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


# Compute and plot the scaling relationship for the surface to volume ratio
phi_range = np.linspace(0.001, 0.08, 100)
delta = 0.024
phi_max = 0.078
wmin = 0.5
alpha = 3.3
Lam = 12 * alpha * delta / (3 * alpha - 1)
dphi_max = phi_range - phi_max

# Compute and plot the scaling relationship for the surface to volume ratio
for i, (g, d) in enumerate(singular_errs[
        singular_errs['quantity'] == 'avg_rho_ratio'].groupby(
        ['interval'], sort=False)):
    lower = Lam * d['lower'].values[0]
    upper = Lam * d['upper'].values[0]
    # w_lower = Lam * ( lower * (1 - phi_range)/phi_range + 1)
    # w_upper = Lam * ( upper * (1 - phi_range)/phi_range + 1)
    # w_[w_lower > 1] = 1
    w_lower = wmin + lower * (dphi_max/phi_max**2) * (dphi_max/phi_max - 1)
    w_upper = wmin + upper * (dphi_max/phi_max**2) * (dphi_max/phi_max - 1)
    ax.fill_between(phi_range, w_upper, w_lower, color=cor['light_black'],
                    alpha=0.4)
#
ax.set_ylim([0.5, 1])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_ylabel('w\naverage width [µm]', fontsize=6)
plt.savefig('../../figures/Fig2_width_prediction_wildtype.pdf')
# %%
(1/0.2 * (wmin * Lam**-1 - 1) + 1)**-1
# %%
