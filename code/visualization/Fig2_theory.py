# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

size_emp = pd.read_csv(
    '../../data/mcmc/size_data_empirical_summaries_wide.csv')
ms_emp = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
per_data = pd.read_csv(
    '../../data/mcmc/ppGpp_perturbation_parameter_summaries.csv')
params = pd.read_csv('../../data/mcmc/theory_parameter_summaries_longform.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
drymass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')
# fig, ax = plt.subplots(1, )

fig, ax = plt.subplots(1, 2, figsize=(4, 2))

ax[0].set_ylim([0, 600])
ax[1].set_ylim([100, 400])

for g, d in prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

for g, d in drymass.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    prot = d[d['quantity'] == 'prot_per_cell']
    rho_dry = d[d['quantity'] == 'drymass_density']
    ax[0].hlines(prot['median_value'], lam['2.5%'], lam['97.5%'],
                 color=cor['primary_blue'], linewidth=1)
    # ax[1].hlines(rho_dry['median_value'], lam['2.5%'],
    #  lam['97.5%'], color=cor['primary_blue'], linewidth=1)
    ax[0].vlines(lam['median_value'], prot['2.5%'], prot['97.5%'],
                 color=cor['primary_blue'], linewidth=1)
    # ax[1].vlines(lam['median_value'], rho_dry['2.5%'],
    #  rho_dry['97.5%'], color=cor['primary_blue'], linewidth=1)
    ax[0].plot(lam['median_value'], prot['median_value'], 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
               markeredgewidth=1, ms=5)
    # ax[1].plot(lam['median_value'], rho_dry['median_value'], 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
    #    markeredgewidth=1, ms=5)

# fig, ax = plt.subplots(2, 2, figsize=(3, 2), sharex=True)
# plt.tight_layout()


# %%
# Inferred constants
fig, ax = plt.subplots(2, 2, figsize=(3, 3), sharex=True)
ax = ax.ravel()
ax[0].set_ylim([0, 0.25])
ax[1].set_ylim([0, 5])
ax[2].set_ylim([0, 30])
ax[3].set_ylim([0, 300])
ax[0].set_ylabel('$\phi_{mem}$\nmembrane allocation', fontsize=6)
ax[1].set_ylabel(r'$\rho_{mem}$'+'\nmembrane density [fg/µm$^2$]', fontsize=6)
ax[2].set_ylabel('$m_{peri}$\nmembrane protein [fg/cell]', fontsize=6)
ax[3].set_ylabel('$\kappa$\ndensity ratio [µm$^{-1}$]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]\n$\lambda$', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]\n$\lambda$', fontsize=6)
# Plot the inferred constants
interval_colors = {'95%': cor['pale_blue'],
                   '75%': cor['light_blue'],
                   '25%': cor['primary_blue'],
                   'median': cor['blue']}

axes = {'phi_mem_mu': ax[0], 'rho_mem_mu': ax[1],
        'm_peri': ax[2], 'kappa': ax[3]}

for g, d in params[params['quantity'].isin(axes.keys())].groupby(['quantity', 'interval'], sort=False):
    _ax = axes[g[0]]
    inter = np.ones(2)
    _ax.fill_between([0, 2.5], d['lower'].values[0] * inter, d['upper'].values[0] * inter,
                     color=interval_colors[g[1]], alpha=0.5)

# Plot the empirical data
for g, d in ms_data.groupby('dataset_name'):
    d = d[d['localization'] == 'membrane']
    fmt = size.viz.style_point(g, alpha=0.5)
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], **fmt, ms=4)

for g, d in ms_emp.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    rho = d[d['quantity'] == 'ms_rho_mem']
    peri = d[d['quantity'] == 'ms_m_peri']
    kappa = d[d['quantity'] == 'ms_kappa']
    for i, _d in enumerate([rho, peri, kappa]):
        ax[i+1].vlines(_d['growth_rate_hr'], _d['2.5%'],
                       _d['97.5%'], linewidth=0.5, color=fmt['color'])
        ax[i+1].plot(_d['growth_rate_hr'], _d['median_value'], **fmt, ms=4)

# Our data
for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    phi = d[d['quantity'] == 'phi_mem']
    rho = d[d['quantity'] == 'rho_mem']
    peri = d[d['quantity'] == 'm_peri']
    for i, _d in enumerate([phi, rho, peri]):
        if (i == 2) & (g == 'LB'):
            continue
        ax[i].hlines(_d['median_value'], lam['2.5%'], lam['97.5%'],
                     linewidth=1, color=cor['primary_blue'])
        ax[i].vlines(lam['median_value'], _d['2.5%'], _d['97.5%'],
                     linewidth=1, color=cor['primary_blue'])
        ax[i].plot(lam['median_value'], _d['median_value'], 'o',
                   markerfacecolor='w', markeredgecolor=cor['primary_blue'],
                   markeredgewidth=1, ms=4)
plt.tight_layout()
plt.savefig('../../figures/Fig2_theory_constants.pdf', bbox_inches='tight')

# %%
posterior = pd.read_csv(
    '../../data/mcmc/theory_parameter_posterior_samples.csv')
fig, ax = plt.subplots(4, 1, figsize=(1.5, 3))
ax = ax.ravel()
for a in ax:
    a.set_yticks([])
    a.set_facecolor('none')
    a.tick_params(direction='out', color=cor['primary_black'], width=2)
    a.spines['bottom'].set_visible(True)
    a.spines['bottom'].set_color(cor['dark_blue'])
    a.spines['bottom'].set_linewidth(0.2)

ax[0].set_xlim([0.1, 0.14])
ax[0].set_xticks([0.1, 0.12, 0.13, 0.14])
ax[0].set_xlabel('membrane allocation\n$\phi_{mem}$', fontsize=6)
ax[1].set_xlim([1, 4])
ax[1].set_xticks([1, 2, 3, 4])
ax[1].set_xlabel('membrane density [fg/µm$^2$]\n' +
                 r'$\rho_{mem}$', fontsize=6)
ax[2].set_xlim([8, 12])
ax[2].set_xticks([8, 9, 10, 11, 12])
ax[2].set_xlabel('periplasmic protein [fg/cell]\n$m_{peri}$', fontsize=6)
ax[3].set_xlim([75, 175])
ax[3].set_xticks([75, 100, 125, 150, 175])
ax[3].set_xlabel('density ratio [µm$^{-1}$]\n$\kappa$', fontsize=6)

for i, p in enumerate(['phi_mem_mu', 'rho_mem_mu', 'm_peri', 'kappa']):
    ax[i].hist(posterior[p].values, bins=100, color=cor['light_blue'], alpha=0.75,
               edgecolor=cor['dark_blue'], linewidth=0.1)

plt.subplots_adjust(hspace=0.9)
plt.savefig('../../figures/Fig2_param_posteriors.pdf', bbox_inches='tight')

# %%
# Plot the theory
theory = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')
fig, ax = plt.subplots(1, 2, figsize=(4, 2))
ax[0].set_xlabel('ribosomal allocation\n$\phi_{Rb}$', fontsize=6)
ax[1].set_xlabel('ribosomal allocation\n$\phi_{Rb}$', fontsize=6)

ax[0].set_ylabel('$S_{A}/V$\nsurface to volume [µm$^{-1}$]', fontsize=6)
ax[1].set_ylabel('$w$\ncell width [µm]', fontsize=6)
for a in ax:
    a.set_xlim([0.08, 0.15])
interval_colors = {'95%': cor['pale_green'],
                   '75%': cor['light_green'],
                   '25%': cor['primary_green'],
                   'median': cor['blue']}
axes = {'SAV_theory': ax[0], 'width_theory': ax[1]}
for g, d in theory.groupby(['quantity', 'interval'], sort=False):
    axes[g[0]].fill_between(d['phiRb'], d['lower'],
                            d['upper'], alpha=0.5, color=interval_colors[g[1]])


for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['median_value'], d['surface_to_volume'], **fmt, ms=4)
    ax[1].plot(d['median_value'], d['width'], **fmt, ms=4)

for g, d in wt_data.groupby('carbon_source'):
    if g != 'glucose':
        continue
    phi = d[d['quantity'] == 'phi_Rb']
    sav = d[d['quantity'] == 'surface_to_volume']
    width = d[d['quantity'] == 'width']
    for i, _d in enumerate([sav, width]):
        ax[i].hlines(_d['median_value'], phi['2.5%'], phi['97.5%'],
                     linewidth=1, color=cor['primary_blue'])
        ax[i].vlines(phi['median_value'], _d['2.5%'], _d['97.5%'],
                     linewidth=1, color=cor['primary_blue'])
        ax[i].plot(phi['median_value'], _d['median_value'], 'o',
                   markerfacecolor='w', markeredgecolor=cor['primary_blue'],
                   markeredgewidth=1, ms=4)

key = {'relA, 0': cor['light_red'],
       'relA, 1': cor['primary_red'],
       'relA, 2': cor['red'],
       'meshI, 0': cor['light_purple'],
       'meshI, 100': cor['primary_purple']}
for g, d in per_data.groupby(['overexpression', 'inducer_conc']):
    idx = f'{g[0]}, {g[1]}'
    c = key[idx]
    phi = d[d['quantity'] == 'phiRb']
    sav = d[d['quantity'] == 'surface_to_volume']
    width = d[d['quantity'] == 'width']
    for i, _d in enumerate([sav, width]):
        ax[i].hlines(_d['median_value'], phi['2.5%'], phi['97.5%'],
                     linewidth=1, color=c)
        ax[i].vlines(phi['median_value'], _d['2.5%'], _d['97.5%'],
                     linewidth=1, color=c)
        ax[i].plot(phi['median_value'], _d['median_value'], 'o',
                   markerfacecolor='w', markeredgecolor=c,
                   markeredgewidth=1, ms=4)


plt.tight_layout()
# plt.savefig('../../figures/Fig2_theory_pred.pdf', bbox_inches='tight')
