# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ms_emp = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
size_emp = pd.read_csv(
    '../../data/mcmc/size_data_empirical_summaries_wide.csv')
fits = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')
theo = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')
wt = pd.read_csv('../../data/mcmc/wildtype_posterior_parameter_summaries.csv')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3.25, 1.65), sharex=True)
ax[0].set_ylim(0, 0.25)
ax[1].set_ylim(0, 0.125)
ax[0].set_xlim([0.05, 0.35])
ax[1].set_xlim([0.05, 0.35])
ax[1].set_yticks([0, 0.05, 0.1])
ax[0].set_yticks([0, 0.1, 0.2])
for a in ax:
    a.set_xlabel(
        'ribosomal allocation\n  $ \phi_{rib} \propto M_{RNA} / M_{prot}^{(tot)}$', fontsize=6)
ax[0].set_ylabel('$\phi_{mem}$\nmembrane allocation', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic allocation', fontsize=6)
for g, d in ms_data.groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    phiRb = d[d['localization'] == 'ribosomal sector']
    phi_mem = d[d['localization'] == 'membrane']
    phi_peri = d[d['localization'] == 'periplasm']
    ax[0].plot(phiRb['mass_frac'], phi_mem['mass_frac'], **fmt, ms=5)
    ax[1].plot(phiRb['mass_frac'], phi_peri['mass_frac'], **fmt, ms=5)

phi_mem_relation = theo[theo['quantity'] == 'phi_mem_pred']
phi_peri_relation = theo[theo['quantity'] == 'phi_peri_pred']
inter_colors = {'95%': 'pale_', '75%': 'light_',
                '25%': 'primary_', 'median': ''}
c = ['blue', 'purple']
for i, p in enumerate([phi_mem_relation, phi_peri_relation]):
    for g, d in p.groupby('interval', sort=False):
        ax[i].fill_between(d['phiRb'], d['lower'], d['upper'], alpha=0.5,
                           color=cor[f'{inter_colors[g]}{c[i]}'], zorder=1000)


for g, d in wt.groupby(['carbon_source']):
    phiRb = d[d['quantity'] == 'phiRb']
    phi_mem = d[d['quantity'] == 'phi_mem']
    phi_peri = d[d['quantity'] == 'phi_peri']
    for i, p in enumerate([phi_mem, phi_peri]):
        if len(p) > 0:
            ax[i].hlines(p['median_value'], phiRb['2.5%'],
                         phiRb['97.5%'], lw=1, color=cor['primary_black'])
            ax[i].vlines(phiRb['median_value'], p['2.5%'],
                         p['97.5%'], lw=1, color=cor['primary_black'])
            ax[i].plot(phiRb['median_value'],  p['median_value'], 'o',
                       ms=5, markerfacecolor='w', markeredgecolor=cor['primary_black'],
                       markeredgewidth=1)

plt.tight_layout()
plt.savefig('../../figures/Fig2.4_param_plots.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
ax.set_xlim([0.11, 0.7])
sav = theo[theo['quantity'] == 'SAV_theory']

for g, d in sav.groupby('interval', sort=False):
    ax.fill_between(d['phiRb'] / 0.4558, d['lower'], d['upper'],
                    alpha=0.75, color=cor[f'{inter_colors[g]}green'])

for g, d in size_emp.groupby(['source']):
    fmt = size.viz.style_point(g)
    ax.plot(d['median_value'] / 0.4558, d['surface_to_volume'], **fmt, ms=5)

for g, d in wt.groupby('carbon_source'):
    phi = wt[wt['quantity'] == 'phiRb']
    sav = wt[wt['quantity'] == 'surface_to_volume']
    ax.hlines(sav['median_value'], phi['2.5%']/0.4558,
              phi['97.5%']/0.4558, lw=1, color=cor['primary_black'])
    ax.vlines(phi['median_value']/0.4558, sav['2.5%'],
              sav['97.5%'], lw=1, color=cor['primary_black'])
    ax.plot(phi['median_value']/0.4558, sav['median_value'], 'o', ms=5, markerfacecolor='w',
            markeredgecolor=cor['primary_black'], markeredgewidth=1)
ax.legend(fontsize=6, handletextpad=0.1)
ax.set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(tot)}$', fontsize=6)
ax.set_ylabel('$S_A/V$ [µm$^{-1}$]\nsurface-to-volume', fontsize=6)
plt.savefig('../../figures/Fig2.4_sav_prediction.pdf', bbox_inches='tight')
