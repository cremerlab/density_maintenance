# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

size_emp = pd.read_csv(
    '../../data/mcmc/size_data_empirical_summaries_wide.csv')
size_emp = size_emp[size_emp['source'] != 'Taheri-Araghi et al. 2015']
ms_emp = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
params = pd.read_csv('../../data/mcmc/theory_parameter_summaries_longform.csv')
posts = pd.read_csv('../../data/mcmc/theory_parameter_posterior_samples.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
drymass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')
pred = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')

# %%
#
fig, ax = plt.subplots(2, 1, figsize=(2, 2), sharex=True)

# Labels
ax[-1].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel(
    r'$\rho_{dry} \approx \rho_{cyt}$ [fg / fL]' + '\ndrymass density', fontsize=6)
ax[1].set_ylabel(r'$\rho_{mem}$  [fg / µm$^2$]' +
                 '\nmembrane density', fontsize=6)

# Limits
ax[0].set_ylim([100, 450])
ax[1].set_ylim([0, 6])
ax[0].set_xlim([0, 2.2])

for g, d in drymass.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)

for g, d in ms_emp[ms_emp['quantity'].isin(['ms_rho_mem'])].groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['light_blue']
    ax[1].vlines(d['growth_rate_hr'], d['2.5%'],
                 d['97.5%'], linewidth=1, color=fmt['color'])
    ax[1].plot(d['growth_rate_hr'], d['median_value'], **fmt)

axes = {'drymass_mu': [ax[0], cor['dark_black']],
        'rho_mem_mu': [ax[1], cor['blue']]}
for g, d in params[params['quantity'].isin(['drymass_mu', 'rho_mem_mu']) &
                   (params['interval'] == 'median')].groupby('quantity'):
    axes[g][0].hlines(d['lower'], 0, 2.5, linewidth=1, color=axes[g][1])
ax[0].legend(fontsize=6)

plt.tight_layout()
plt.savefig('../../figures/Fig2.1_theory_constants.pdf')
# %%
fig, ax = plt.subplots(2, 1, figsize=(1, 2))
for a in ax:
    a.set_yticks([])
    a.set_facecolor('#FFF')
    a.spines['bottom'].set_visible(True)
    a.spines['bottom'].set_color(cor['dark_black'])
    a.spines['bottom'].set_linewidth(0.2)
ax[0].set_xlim([275, 300])
ax[0].set_xticks([280, 290, 300])
ax[1].set_xlim([2, 3.5])
ax[1].set_xticks([2, 2.5, 3, 3.5])
ax[0].set_xlabel('drymass density\n' +
                 r'$\rho_{dry}\approx\rho_{cyt}$ [fg / fL]', fontsize=6)
ax[1].set_xlabel('membrane density\n' + r'$\rho_{mem}$ [fg / µm]', fontsize=6)
_ = ax[0].hist(posts['drymass_mu'], bins=50,
               color=cor['light_black'], edgecolor=cor['primary_black'])
_ = ax[1].hist(posts['rho_mem_mu'], bins=50,
               color=cor['light_blue'], edgecolor=cor['primary_blue'])
plt.tight_layout()
plt.savefig('../../figures/Fig2.1_theory_constants_dists.pdf')


# %%
fig, ax = plt.subplots(1, 1, figsize=(1, 1))
ax.set_facecolor('#FFF')
ax.spines['bottom'].set_visible(True)
ax.spines['bottom'].set_linewidth(0.2)
ax.spines['bottom'].set_color(cor['dark_black'])
ax.set_yticks([])
ax.set_xlim([75, 150])
ax.set_xlabel('density ratio\n$\kappa$ [µm$^{-1}$]', fontsize=6)
_ = ax.hist(posts['kappa'], bins=50,
            color=cor['light_red'], edgecolor=cor['red'])
# plt.tight_layout()
plt.savefig('../../figures/Fig2.1_kappa_dist.pdf')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)
for a in ax:
    a.set_xlabel('ribosomal protien\nallocation $\phi_{Rb}$', fontsize=6)
ax[0].set_ylabel('$S_A/V$ [µm$^{-1}$]\naverage surface-to-volume', fontsize=6)
ax[1].set_ylabel('$w$[µm]\naverage cell width', fontsize=6)
ax[0].set_xlim([0.05, 0.3])
ax[0].set_ylim([3, 10])
ax[1].set_ylim([0.4, 1.2])

axes = {'SAV_theory': ax[0], 'width_theory': ax[1]}
inter_colors = {'95%': cor['pale_green'], '75%': cor['light_green'],
                '25%': cor['primary_green'], 'median': cor['green']}
for g, d in pred.groupby(['quantity', 'interval'], sort=False):
    _ax = axes[g[0]].fill_between(
        d['phiRb'], d['lower'], d['upper'], color=inter_colors[g[1]])

for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].hlines(d['surface_to_volume'], d['2.5%'],
                 d['97.5%'], linewidth=1, color=fmt['color'])
    ax[1].hlines(d['width'], d['2.5%'], d['97.5%'],
                 linewidth=1, color=fmt['color'])
    ax[0].plot(d['median_value'], d['surface_to_volume'], **fmt)
    ax[1].plot(d['median_value'], d['width'], **fmt)


for g, d in wt_data.groupby('carbon_source'):
    phi = d[d['quantity'] == 'phiRb']
    width = d[d['quantity'] == 'width']
    sav = d[d['quantity'] == 'surface_to_volume']
    for i, q in enumerate([sav, width]):
        ax[i].vlines(phi['median_value'], q['2.5%'], q['97.5%'],
                     linewidth=1, color=cor['primary_blue'], label='__nolegend__')
        ax[i].hlines(q['median_value'], phi['2.5%'], phi['97.5%'],
                     linewidth=1, color=cor['primary_blue'], label='__nolegend__')
        ax[i].plot(phi['median_value'], q['median_value'], 'o', markerfacecolor='w',
                   markeredgecolor=cor['primary_blue'], markeredgewidth=1, ms=5, label='__nolegend__')

ax[0].plot([], [], 'o', markerfacecolor='w', color=cor['primary_blue'],
           markeredgewidth=1, markeredgecolor=cor['primary_blue'], label='This study')
ax[0].legend(fontsize=6)
plt.savefig('../../figures/Fig2.1_pred_data.pdf')
# plt.tight_layout()
