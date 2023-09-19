# %%
import numpy as np
import pandas as pd
import size.viz
import matplotlib.pyplot as plt
cor, pal = size.viz.matplotlib_style()

# Load the necessary datasets
ms_emp = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
pred = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')
pars = pd.read_csv('../../data/mcmc/theory_parameter_summaries_wide.csv')
posts = pd.read_csv('../../data/mcmc/theory_parameter_posterior_samples.csv')

fig, ax = plt.subplots(2, 1, figsize=(2, 2), sharex=True)
ax[0].set_ylim([0, 30])
ax[0].set_ylabel('periplasmic protein\nmass [fg/cell]\n', fontsize=6)
ax[1].set_ylabel('total protein\nmass [fg/cell]', fontsize=6)
ax[1].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
# ax[1].set_ylim([80, 1000])
# ax[1].set_yscale('log')
ax[1].set_yticks([100, 250, 500, 1000])
ax[1].set_yticklabels([100, 250, 500, 1000])

_pars = pars[pars['quantity'] == 'm_peri']
ax[0].fill_between([0, 2.5], _pars['2.5%'].values * [1, 1], _pars['97.5%'].values
                   * [1, 1], color=cor['pale_purple'])
ax[0].hlines(_pars['median_value'], 0, 2.5, lw=1,
             color=cor['purple'])


# Plot the literature data
for g, d in ms_emp[ms_emp['quantity'] == 'ms_m_peri'].groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    fmt['color'] = cor['light_purple']
    ax[0].vlines(d['growth_rate_hr'], d['2.5%'], d['97.5%'],
                 linewidth=1, color=fmt['color'])
    ax[0].plot(d['growth_rate_hr'], d['median_value'], **fmt)

for g, d in prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

int_color = {'95%': cor['pale_black'], 'median': cor['black']}
for g, d in pred[(pred['quantity'] == 'pred_lam_prot') &
                 (pred['interval'].isin(['median', '95%']))].groupby('interval', sort=False):
    ax[1].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                       color=int_color[g])

# Plot the inferred parameter sumamries
for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    m_peri = d[d['quantity'] == 'm_peri']
    prot = d[d['quantity'] == 'prot_per_cell']

    for i, _d in enumerate([m_peri, prot]):
        if len(_d) == 0:
            continue
        ax[i].hlines(_d['median_value'], lam['2.5%'],
                     lam['97.5%'], lw=1, color=cor['primary_blue'])
        ax[i].vlines(lam['median_value'], _d['2.5%'],
                     _d['97.5%'], lw=1, color=cor['primary_blue'])
        ax[i].plot(lam['median_value'], _d['median_value'], 'o', markerfacecolor='w',
                   markeredgecolor=cor['primary_blue'], markeredgewidth=1, ms=4)
plt.savefig('../../figures/Figx_periplasm_protein_trends.pdf')

# %%
fig = plt.figure(figsize=(2, 2))
gs = fig.add_gridspec(2, 2)
ax0 = fig.add_subplot(gs[0, :])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[1, 1])
ax = [ax0, ax1, ax2]
ax[0].set_xlim([6, 14])
ax[1].set_xlim([75, 150])
ax[2].set_xlim([1.5, 3])

for a in ax:
    a.set_facecolor('none')
    a.spines['bottom'].set_visible(True)
    a.spines['bottom'].set_linewidth(0.2)
    a.spines['bottom'].set_color(cor['primary_black'])
    a.set_yticks([])


ax[0].set_xlabel(
    'periplasmic protein\n$M_{prot}^{(peri)}$ [fg / cell]', fontsize=6)
ax[1].set_xlabel(
    'minimum protein mass\n$M_{prot_0}^{(tot)}}$ [fg / cell]', fontsize=6)
ax[2].set_xlabel('protein mass constant\n$k$ [hr]', fontsize=6)


_ = ax[0].hist(posts['m_peri'], bins=100,
               color=cor['pale_purple'], edgecolor=cor['primary_purple'])
_ = ax[1].hist(np.exp(posts['log_prot_intercept']), bins=100,
               color=cor['pale_black'], edgecolor=cor['primary_black'])
_ = ax[2].hist(np.exp(posts['log_prot_slope']), bins=100,
               color=cor['pale_black'], edgecolor=cor['primary_black'])
plt.subplots_adjust(hspace=0.5)
plt.savefig('../../figures/FigX_periplasm_protein_posteriors.pdf',
            bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.75), sharex=True)
ax[0].set_xlim([0, 1.5])
for a in ax:
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\phi_{peri}$\nperiplasmic allocation', fontsize=6)
ax[1].set_ylabel(r'$\rho_{peri}$ [fg/fL]' +
                 '\nperiplasmic allocation', fontsize=6)


interval_colors = {'95%': cor['pale_green'], '75%': cor['light_green'],
                   '25%': cor['primary_green'], 'median': cor['green']}
axes = {'phi_peri_pred': ax[0], 'rho_peri_pred': ax[1]}
for g, d in pred[pred['quantity'].isin(['phi_peri_pred', 'rho_peri_pred'])].groupby(['quantity', 'interval'], sort=False):
    axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                            color=interval_colors[g[1]], alpha=0.75)

for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], **fmt)

for g, d in ms_emp.groupby('source'):
    _d = d[d['quantity'] == 'ms_rho_peri']
    fmt = size.viz.style_point(g)
    ax[1].vlines(_d['growth_rate_hr'], _d['2.5%'],
                 _d['97.5%'], lw=1, color=fmt['color'])
    ax[1].plot(_d['growth_rate_hr'], _d['median_value'], **fmt)


for g, d in wt_data.groupby('carbon_source'):
    if g == 'LB':
        continue
    lam = d[d['quantity'] == 'growth_rate']
    phi_peri = d[d['quantity'] == 'phi_peri']
    rho_peri = d[d['quantity'] == 'rho_peri']
    for i, _d in enumerate([phi_peri, rho_peri]):
        ax[i].hlines(_d['median_value'], lam['2.5%'],
                     lam['97.5%'], lw=1, color=cor['primary_blue'])
        ax[i].vlines(lam['median_value'], _d['2.5%'],
                     _d['97.5%'], lw=1, color=cor['primary_blue'])
        ax[i].plot(lam['median_value'], _d['median_value'], 'o', markerfacecolor='w',
                   markeredgecolor=cor['primary_blue'], markeredgewidth=1, ms=4)
plt.tight_layout()
plt.savefig('../../figures/FigX_periplasm_allocation_density_prediction.pdf',
            bbox_inches='tight')
