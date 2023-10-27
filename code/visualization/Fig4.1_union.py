# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

vol_data = pd.read_csv(
    '../../data/literature/collated_literature_volume_data.csv')
vol_data = vol_data[~vol_data['source'].isin(
    ['Si et al. 2017', 'Taheri-Araghi et al. 2015', 'Basan et al. 2015'])]
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Taheri-Araghi et al. 2015']
pred = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')
pred = pred[pred['quantity'].isin(['SAV_theory', 'width_theory'])]
lam_pred = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

size_emp = pd.read_csv(
    '../../data/mcmc/size_data_empirical_summaries_wide.csv')
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
params = pd.read_csv('../../data/mcmc/theory_parameter_summaries_wide.csv')

# %%
fig, ax = plt.subplots(2, 1, figsize=(2, 2.5))
ax[0].set_xlim([0.1, 0.7])
ax[0].set_ylim([0.45, 1.1])
ax[1].set_xlim([0, 2])
ax[1].set_ylim([-0.1, 3.5])
ax[0].set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(cyt)}$', fontsize=6)
ax[0].set_ylabel('w [µm]\ncell width', fontsize=6)
ax[1].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[1].set_ylabel('V [fL]\ncell width', fontsize=6)

width_theo = pred[pred['quantity'] == 'width_theory']
inter_colors = {'95%': cor['pale_green'], '75%': cor['light_green'],
                '25%': cor['primary_green'], 'median': cor['green']}

for g, d in width_theo.groupby('interval', sort=False):
    ax[0].fill_between(d['phiRb']/0.4558, d['lower'], d['upper'], alpha=0.75,
                       color=inter_colors[g])

for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['median_value']/0.4558, d['width'], **fmt)

for g, d in vol_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    phi = d[d['quantity'] == 'phiRb']
    width = d[d['quantity'] == 'width']
    vol = d[d['quantity'] == 'volume']
    lam = d[d['quantity'] == 'growth_rate']
    ax[0].vlines(phi['median_value']/0.4558, width['2.5%'],
                 width['97.5%'], lw=1, color=cor['primary_black'])
    ax[0].hlines(width['median_value'], phi['2.5%']/0.4558,
                 phi['97.5%']/0.4558, lw=1, color=cor['primary_black'])
    ax[0].plot(phi['median_value']/0.4558, width['median_value'], 'o', markeredgewidth=1,
               markeredgecolor=cor['primary_black'], color='w', ms=6)
    ax[1].vlines(lam['median_value'], vol['2.5%'],
                 vol['97.5%'], lw=1, color=cor['primary_black'])
    ax[1].hlines(vol['median_value'], lam['2.5%'],
                 lam['97.5%'], lw=1, color=cor['primary_black'])
    ax[1].plot(lam['median_value'], vol['median_value'], 'o', markeredgewidth=1,
               markeredgecolor=cor['primary_black'], color='w', ms=6)

vol = lam_pred[lam_pred['quantity'] == 'volume_pred']
for g, d in vol.groupby('interval', sort=False):
    ax[1].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], alpha=0.75,
                       color=inter_colors[g])
plt.savefig('../../figures/Fig4.1_width_pred.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_ylim([1, 6])
ax.set_xlim([0.1, 0.7])
ax.set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(cyt)}$', fontsize=6)
ax.set_ylabel(r'$\alpha$'+'\naspect ratio', fontsize=6)
for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    ell = size_data[size_data['source'] == g]
    ax.plot(d['median_value']/0.4558, ell['length_um']/ell['width_um'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    phi = d[d['quantity'] == 'phiRb']
    alpha = d[d['quantity'] == 'aspect_ratio']
    ax.vlines(phi['median_value']/0.4558, alpha['2.5%'],
              alpha['97.5%'], lw=1, color=cor['primary_black'])
    ax.hlines(alpha['median_value'], phi['2.5%']/0.4558,
              phi['97.5%']/0.4558, lw=1, color=cor['primary_black'])
    ax.plot(phi['median_value']/0.4558, alpha['median_value'], 'o', markeredgewidth=1,
            markeredgecolor=cor['primary_black'], color='w', ms=6)

alpha = params[params['quantity'] == 'alpha_mu']
ax.plot([0.05, 0.7], np.ones(2) * alpha['median_value'].values[0],
        color=cor['dark_red'], zorder=1001)
ax.fill_between([0.05, 0.7], alpha['2.5%'], alpha['97.5%'],
                color=cor['pale_red'], alpha=0.5, zorder=1000)
plt.savefig('../../figures/Fig4.1_alpha.pdf')


# %%
fig, ax = plt.subplots(1, 3, figsize=(4.5, 1.75), sharex=True)
ax[0].set_xlim([-0.1, 2.2])
ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 4.5])
ax[2].set_ylim([0, 3])
ax[0].set_ylabel('w [µm]\nwidth', fontsize=6)
ax[1].set_ylabel('$\ell$ [µm]\n length', fontsize=6)
ax[2].set_ylabel('$V$ [µm$^{3}$]\n volume', fontsize=6)
for a in ax:
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    for i, k in enumerate(['width_um', 'length_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[k], **fmt, ms=4)

for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    w = d[d['quantity'] == 'width']
    ell = d[d['quantity'] == 'length']
    vol = d[d['quantity'] == 'volume']
    for i, p in enumerate([w, ell, vol]):
        ax[i].vlines(lam['median_value'], p['2.5%'],
                     p['97.5%'], lw=1, color=cor['primary_black'])
        ax[i].hlines(p['median_value'], lam['2.5%'], lam['97.5%'], lw=1,
                     color=cor['primary_black'])
        ax[i].plot(lam['median_value'], p['median_value'], 'o', markeredgewidth=1,
                   markeredgecolor=cor['primary_black'], color='w', ms=4)

axes = {'width_pred': ax[0], 'length_pred': ax[1], 'volume_pred': ax[2]}
for g, d in lam_pred[lam_pred['quantity'].isin(list(axes.keys()))].groupby(['quantity', 'interval'], sort=False):
    _ax = axes[g[0]]
    _ax.fill_between(d['growth_rate_hr'], d['lower'],
                     d['upper'], alpha=0.75, color=inter_colors[g[1]])
plt.tight_layout()
plt.savefig('../../figures/Fig4.1_lam_preds.pdf')
