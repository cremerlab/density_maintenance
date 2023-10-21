# %%
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
pred = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')
pred = pred[pred['quantity'].isin(['SAV_theory', 'width_theory'])]
emp = pd.read_csv('../../data/mcmc/size_data_empirical_summaries_wide.csv')
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
posts = pd.read_csv('../../data/mcmc/ppGpp_perturbation_samples.csv')
params = pd.read_csv(
    '../../data/mcmc/ppGpp_perturbation_parameter_summaries.csv')

# %%
fig, ax = plt.subplots(5, 4, figsize=(3.5, 1.25))
for a in ax.ravel():
    a.set_yticks([])
    a.set_facecolor('none')
    a.grid(False)
    a.spines['bottom'].set_visible(True)
    a.spines['bottom'].set_linewidth(0.2)
    a.spines['bottom'].set_color(cor['primary_black'])
    a.set_xticks([])
for i in range(5):
    ax[i, 0].set_xlim([0.5, 1.1])
    ax[i, 1].set_xlim([0.2, 0.35])
    ax[i, 2].set_xlim([5.5, 8.5])
    ax[i, 3].set_xlim([0.4, 0.9])


ax[-1, 0].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[-1, 1].set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(tot)}$', fontsize=6)
ax[-1, 2].set_xlabel('surface-to-volume\n$S_A/V$ [µm$^{-1}$]', fontsize=6)
ax[-1, 3].set_xlabel('width\n$w$ [µm]', fontsize=6)
ax[-1, 0].set_xticks([0.5, 0.75, 1])
ax[-1, 1].set_xticks([0.2, 0.25, 0.3, 0.35])
ax[-1, 2].set_xticks([5.5, 7, 8.5])
ax[-1, 3].set_xticks([0.4, 0.6, 0.8])


axes = {'meshI': {100: [0, cor['dark_red'], cor['dark_black']], 0: [1, cor['red'], cor['dark_black']]},
        'relA': {0: [2, cor['primary_red'], cor['dark_black']], 1: [3, cor['light_red'], cor['dark_black']], 2: [4, cor['pale_red'], cor['dark_black']]}}
pars = {'growth_rate': 0, 'phiRb': 1, 'surface_to_volume': 2, 'width': 3}
for g, d in posts[posts['quantity'].isin(list(pars.keys()))].groupby(['quantity', 'overexpression', 'inducer_conc']):
    if g[0] == 'phiRb':
        div = 0.4558
    else:
        div = 1
    _ax = ax[axes[g[1]][g[2]][0], pars[g[0]]]
    _ = _ax.hist(d['value'] / div, bins=100, color=axes[g[1]][g[2]][1],
                 edgecolor=axes[g[1]][g[2]][2], linewidth=0.2, histtype='stepfilled',
                 alpha=0.75)


plt.subplots_adjust(hspace=-0.25)
plt.savefig('../../figures/Fig3.2_perturbation_joyplot.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(1.75, 1))
ax.set_ylim([5.5, 8.5])
ax.set_xlim([0.19, 0.325])
ax.set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(tot)}$', fontsize=6)
ax.set_ylabel('$S_A/V$ [µm$^{-1}$]\nsurface-to-volume', fontsize=6)

# ax[1].set_xlim([0.2, 0.3])

colors = {'meshI': {'m': 's', 100: 'pale_', 0: 'light_'},
          'relA': {'m': 'D', 0: 'primary_', 1: '', 2: 'dark_'}}
for g, d in params.groupby(['overexpression', 'inducer_conc']):
    sav = d[d['quantity'] == 'surface_to_volume']
    phi = d[d['quantity'] == 'phiRb']
    m = colors[g[0]]['m']
    sav_c = cor[f'{colors[g[0]][g[1]]}red']
    ax.hlines(sav['median_value'], phi['2.5%']/0.4558,
              phi['97.5%']/0.4558, lw=1, color=sav_c)
    ax.vlines(phi['median_value']/0.4558, sav['2.5%'],
              sav['97.5%'], lw=1, color=sav_c)
    ax.plot(phi['median_value']/0.4558, sav['median_value'], marker=m, markerfacecolor='w', markeredgecolor=sav_c,
            markeredgewidth=1, ms=4)

plt.savefig('../../figures/Fig3.2_sav_points.pdf', bbox_inches='tight')
