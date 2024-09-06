# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

densities = pd.read_csv(
    '../../data/mcmc/inner_outer_membrane_densities_wide.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')


fig, ax = plt.subplots(2, 2, figsize=(3, 2))
ax[0, 0].set_ylim([0, 0.15])
ax[1, 0].set_ylim([0, 0.15])
ax[0, 1].set_ylim([0, 7])
ax[1, 1].set_ylim([0, 7])
for g, d in ms_data[ms_data['localization'] == 'inner membrane'].groupby('dataset_name'):
    fmt = size.viz.style_point(g, alpha=0.45)
    fmt['color'] = cor['pale_blue']
    ax[0, 0].plot(d['growth_rate_hr'], d['mass_frac'], **fmt, ms=4)

for g, d in ms_data[ms_data['localization'] == 'outer membrane'].groupby('dataset_name'):
    fmt = size.viz.style_point(g, alpha=0.45)
    fmt['color'] = cor['blue']
    ax[1, 0].plot(d['growth_rate_hr'], d['mass_frac'], **fmt, ms=4)

axes = {'ms_rho_mem_inner': [ax[0, 1], cor['pale_blue']],
        'ms_rho_mem_outer': [ax[1, 1], cor['dark_blue']]}
for g, d in densities.groupby(['quantity', 'source']):
    fmt = size.viz.style_point(g[1], alpha=0.45)
    _axes = axes[g[0]]
    fmt['color'] = _axes[1]
    _axes[0].vlines(d['growth_rate_hr'], d['2.5%'],
                    d['97.5%'], lw=0.5, color=fmt['color'], alpha=0.45)
    _axes[0].plot(d['growth_rate_hr'], d['median_value'], **fmt, ms=4)

for a in ax.ravel():
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[0, 0].set_ylabel(
    '$\phi_{mem}^{(inner)}$\ninner membrane\nprotein fraction', fontsize=6)
ax[1, 0].set_ylabel(
    '$\phi_{mem}^{(outer)}$\nouter membrane\nprotein fraction', fontsize=6)
ax[0, 1].set_ylabel(r'$\rho_{mem}^{(inner)}$ [fg / µm$^2$]' +
                    '\ninner membrane\nprotein density', fontsize=6)
ax[1, 1].set_ylabel(r'$\rho_{mem}^{(outer)}$ [fg / µm$^2$]' +
                    '\nouter membrane\nprotein density', fontsize=6)
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(-0.05, 0.95, '(A)', fontsize=6)
fig.text(0.45, 0.95, '(B)', fontsize=6)
fig.text(-0.05, 0.45, '(C)', fontsize=6)
fig.text(0.45, 0.45, '(D)', fontsize=6)
plt.savefig('../../FigS3_inner_outer_membrane_densities.pdf',)
