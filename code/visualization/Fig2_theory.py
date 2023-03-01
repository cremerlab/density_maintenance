# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
delta = 0.029
# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
# Restrict to wildtype
params = params[(params['strain'] == 'wildtype') &
                (params['overexpression'] == 'none') &
                (params['inducer_conc'] == 0)]

singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']

# %%
fig, ax = plt.subplots(4, 1, figsize=(2.25, 3), sharex=True)

ycoords = {'glucoseCAA': 5, 'glucose': 4,
           'glycerol': 3, 'sorbitol': 2, 'acetate': 1}
axes = {'alpha': 0, 'rho_cyt': 1, 'rho_peri': 2, 'rho_ratio': 3}

err_widths = {'95%': 0.5, '75%': 1.5, '25%': 2.5}
for g, d in errs[errs['quantity'].isin(list(axes.keys()))
                 ].groupby(['carbon_source', 'quantity', 'interval']):
    if g[1] in ['rho_cyt', 'rho_peri']:
        mult = 1E6
    else:
        mult = 1
    ax[axes[g[1]]].vlines(ycoords[g[0]], mult * d['lower'], mult * d['upper'], lw=err_widths[g[-1]],
                          color=cor['primary_blue'])


for g, d in medians[medians['quantity'].isin(list(axes.keys()))
                    ].groupby(['carbon_source', 'quantity']):
    if g[1] in ['rho_cyt', 'rho_peri']:
        mult = 1E6
    else:
        mult = 1
    ax[axes[g[1]]].plot(ycoords[g[0]], mult * d['lower'], 'o', markeredgewidth=0.5, ms=3,
                        markeredgecolor=cor['primary_blue'], markerfacecolor='white')

# Add labels
ax[-1].set_xticks(list(ycoords.values()))
ax[0].set_xlim([0.5, 5.5])
labels = ['glucose\n+ CAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']

# Update limits
ax[-1].set_xticklabels(labels)
ax[0].set_ylim([2, 5])
ax[0].set_yticks([2, 3, 4, 5])
ax[1].set_ylim([0.4, 0.8])
ax[1].set_yticks([0.5, 0.6, 0.7])
ax[2].set_ylim([0, 0.2])
ax[2].set_yticks([0, 0.1, 0.2])
ax[3].set_ylim([0, 0.3])
ax[3].set_yticks([0, 0.1, 0.2, 0.3])

# Add labels
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)
plt.savefig('../../figures/Fig2_constant_parameters.pdf', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))

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
phi_range = np.linspace(0.001, 0.06, 100)
phi_max = 0.1
wmin = 0.5
alpha = 3.1
Lam = 12 * alpha * delta / (3 * alpha - 1)
dphi_max = phi_range - phi_max

# Compute and plot the scaling relationship for the surface to volume ratio
for i, (g, d) in enumerate(singular_errs[
        singular_errs['quantity'] == 'avg_rho_ratio'].groupby(
        ['interval'], sort=False)):
    lower = 12 * 3 * delta * d['lower'].values[0] / (3 * 3 - 1)
    upper = 12 * 3 * delta * d['upper'].values[0] / (3 * 3 - 1)
    w_lower = wmin + lower * (dphi_max/phi_max**2) * (dphi_max/phi_max - 1)
    w_upper = wmin + upper * (dphi_max/phi_max**2) * (dphi_max/phi_max - 1)
    ax.fill_between(phi_range, w_upper, w_lower, color=cor['light_black'],
                    alpha=0.4)

ax.set_ylim([0.5, 1])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_ylabel('w\naverage width [µm]', fontsize=6)
plt.savefig('../../figures/Fig2_width_prediction_wildtype.pdf')
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

for g, d in errs[errs['quantity'] == 'rho_ratio'].groupby(['carbon_source', 'interval']):

    lam = growth_rates[growth_rates['carbon_source']
                       == g[0]]['growth_rate_hr'].values[0]
    if g[0] == 'LB':
        continue
    ax.vlines(lam, d['lower'], d['upper'], linewidth=err_widths[g[1]],
              color=cor['primary_blue'])

for g, d in medians.groupby(['carbon_source']):
    if g == 'LB':
        continue
    lam = growth_rates[growth_rates['carbon_source']
                       == g]['growth_rate_hr'].values[0]
    ax.plot(lam, d[d['quantity'] == 'rho_ratio']['lower'], 'D', ms=4,
            markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
            markerfacecolor='white')
