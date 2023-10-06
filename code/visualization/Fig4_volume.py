# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Taheri-Araghi et al. 2015']
vol_data = pd.read_csv(
    '../../data/literature/collated_literature_volume_data.csv')
vol_data = vol_data[vol_data['source'] != 'Taheri-Araghi et al. 2015']
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
pert_data = pd.read_csv(
    '../../data/mcmc/ppGpp_perturbation_parameter_summaries.csv')
pred = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')
posts = pd.read_csv('../../data/mcmc/theory_parameter_posterior_samples.csv')
pars = pd.read_csv('../../data/mcmc/theory_parameter_summaries_wide.csv')

# %%
fig, ax = plt.subplots(2, 2, figsize=(3.5, 3))
ax = ax.ravel()
ax[0].set_ylim([0, 6])
ax[2].set_ylim([1, 4.5])
ax[3].set_ylim([-0.2, 4])
ax[0].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[2].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[1].set_xlabel('aspect ratio\n'+r'$\alpha$', fontsize=6)
ax[0].set_ylabel('$\ell/w$\naspect ratio', fontsize=6)
ax[2].set_ylabel('$\ell$ [µm]\naverage length', fontsize=6)
ax[3].set_ylabel('$V$ [µm$^3$]\naverage volume', fontsize=6)
ax[0].set_xlim([0, 2.1])
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['length_um'] / d['width_um'], **fmt)
    ax[2].plot(d['growth_rate_hr'], d['length_um'], **fmt)
    ax[3].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)

for g, d in vol_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[3].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    alpha = d[d['quantity'] == 'aspect_ratio']
    ell = d[d['quantity'] == 'length']
    vol = d[d['quantity'] == 'volume']
    for i, q in enumerate([alpha, ell, vol]):
        if i > 0:
            i += 1
        ax[i].hlines(q['median_value'], lam['2.5%'], lam['97.5%'],
                     lw=1, color=cor['primary_blue'])
        ax[i].vlines(lam['median_value'], q['2.5%'], q['97.5%'],
                     lw=1, color=cor['primary_blue'])
        ax[i].plot(lam['median_value'], q['median_value'], 'o', markerfacecolor='w',
                   markeredgecolor=cor['primary_blue'], markeredgewidth=1, ms=4)

alpha = pars[pars['quantity'] == 'alpha_mu']
ax[0].fill_between([0, 2.6], alpha['2.5%'].values * [1, 1], alpha['97.5%'].values * [1, 1],
                   color=cor['pale_black'])
ax[0].hlines(alpha['median_value'], 0, 2.6,
             color=cor['primary_black'], linewidth=1)
axes = {'length_pred': ax[2], 'volume_pred': ax[3]}
cred_colors = {'95%': cor['pale_green'], '75%': cor['light_green'],
               '25%': cor['primary_green'], 'median': cor['green']}
for g, d in pred[pred['quantity'].isin(['length_pred', 'volume_pred'])].groupby(['quantity', 'interval'], sort=False):
    _ax = axes[g[0]]
    _ax.fill_between(d['growth_rate_hr'], d['lower'],
                     d['upper'], color=cred_colors[g[1]], alpha=0.75)

_ = ax[1].hist(posts['alpha_mu'], bins=100,
               color=cor['light_black'], edgecolor=cor['black'])
ax[1].set_facecolor('none')
ax[1].set_yticks([])
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['bottom'].set_linewidth(0.2)
ax[1].spines['bottom'].set_color(cor['primary_black'])
# plt.tight_layout()
plt.savefig('../../figures/Fig4_aspect_ratio_volume_pred.pdf')

# %%
