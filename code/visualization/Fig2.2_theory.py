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

fig, ax = plt.subplots(2, 2,  figsize=(3.5, 2))
ax = ax.ravel()
ax[1].set_ylim([0, 0.35])
ax[2].set_ylim([0, 0.30])
ax[3].set_ylim([0, 0.13])
ax[3].set_yticks([0, 0.05, 0.10])
ax[0].axis(False)

for a in ax[1:]:
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[1].set_ylabel('$\phi_{Rb}$\nribosomal allocation', fontsize=6)
ax[2].set_ylabel('$\phi_{mem}$\nmembrane allocation', fontsize=6)
ax[3].set_ylabel('$\phi_{peri}$\nperiplasmic allocation', fontsize=6)

inter_colors = {'95%': 'pale_', '75%': 'light_',
                '25%': 'primary_', 'median': ''}

# Plot the mass spectrometry measurements
for g, d in ms_data.groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax[0].plot([], [], **fmt, ms=4)
    phiRb = d[d['localization'] == 'ribosomal sector']
    phi_mem = d[d['localization'] == 'membrane']
    phi_peri = d[d['localization'] == 'periplasm']
    for i, _d in enumerate([phiRb, phi_mem, phi_peri]):
        ax[i+1].plot(_d['growth_rate_hr'], _d['mass_frac'], **fmt, ms=4)

# Plot the parameter fits
for g, d in fits.groupby('interval', sort=False):
    phi_Rb = d[d['quantity'] == 'phi_Rb_pred']
    phi_mem = d[d['quantity'] == 'phi_mem_pred']
    phi_peri = d[d['quantity'] == 'phi_peri_pred']
    for i, (_d, _c) in enumerate(zip([phi_Rb, phi_mem, phi_peri], ['gold', 'blue', 'purple'])):
        ax[i+1].fill_between(_d['growth_rate_hr'], _d['lower'], _d['upper'],
                             color=cor[f'{inter_colors[g]}{_c}'], alpha=0.75,
                             zorder=1000)

# Plot our measurements
for g, d in wt.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    phi_Rb = d[d['quantity'] == 'phiRb']
    phi_mem = d[d['quantity'] == 'phi_mem']
    phi_peri = d[d['quantity'] == 'phi_peri']
    for i, _d in enumerate([phi_Rb, phi_mem, phi_peri]):
        if len(_d) == 0:
            continue
        ax[i+1].hlines(_d['median_value'], lam['2.5%'],
                       lam['97.5%'], lw=0.75, color=cor['primary_black'])
        ax[i+1].vlines(lam['median_value'], _d['2.5%'],
                       _d['97.5%'], lw=0.75, color=cor['primary_black'])
        ax[i+1].plot(lam['median_value'], _d['median_value'], 'o', markeredgewidth=0.75,
                     markerfacecolor='w', markeredgecolor=cor['primary_black'], ms=4)
ax[0].legend(fontsize=6)
plt.savefig('../../figures/Fig2.2_model_parameters.pdf', bbox_inches='tight')

# %%
# Plot the predictions
fig, ax = plt.subplots(2, 1, figsize=(2, 3.5))
ax[0].set_ylim([3.5, 8.5])
ax[1].set_ylim([0.5, 1.3])
ax[0].set_ylabel('$S_A/V$ [µm$^{-1}$]\naverage surface-to-volume', fontsize=6)
ax[1].set_ylabel('w [µm]\naverage cell width', fontsize=6)
for a in ax:
    a.set_xlim([0.05, 0.35])
    a.set_xlabel('ribosomal allocation\n$\phi_{Rb}$', fontsize=6)

axes = {'SAV_theory': ax[0], 'width_theory': ax[1]}
for g, d in theo.groupby(['interval', 'quantity'], sort=False):
    axes[g[1]].fill_between(d['phiRb'], d['lower'], d['upper'],
                            color=cor[f'{inter_colors[g[0]]}green'], alpha=0.75)

for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['median_value'], d['surface_to_volume'], **fmt, ms=5)
    ax[1].plot(d['median_value'], d['width'], **fmt, ms=5)

for g, d in wt.groupby('carbon_source'):
    phi_Rb = d[d['quantity'] == 'phiRb']
    sav = d[d['quantity'] == 'surface_to_volume']
    width = d[d['quantity'] == 'width']
    for i, _d in enumerate([sav, width]):
        ax[i].vlines(phi_Rb['median_value'], _d['2.5%'], _d['97.5%'], linewidth=1,
                     color=cor['primary_black'])
        ax[i].hlines(_d['median_value'], phi_Rb['2.5%'], phi_Rb['97.5%'], linewidth=1,
                     color=cor['primary_black'])
        ax[i].plot(phi_Rb['median_value'], _d['median_value'], 'o', ms=5,
                   markerfacecolor='w', markeredgecolor=cor['primary_black'], markeredgewidth=1)
ax[0].legend(fontsize=6)
plt.savefig('../../figures/Fig2.2_model_predictions.pdf', bbox_inches='tight')
