# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()

phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')

wt_data = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')

theo = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')

# %%
w_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_data['width'] = w_popt[1] + w_popt[0] * phiRb_data['growth_rate_hr']

fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
for g, d in phiRb_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['mass_fraction'] / 0.4558, d['width'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    phi = d[d['quantity'] == 'phi_Rb']
    w = d[d['quantity'] == 'width']

    ax.vlines(phi['median_value'].values / 0.4558, w['2.5%'],
              w['97.5%'], linewidth=1, color=cor['primary_blue'])
    ax.hlines(w['median_value'].values, phi['2.5%'] / 0.4558,
              phi['97.5%'] / 0.4558, linewidth=1, color=cor['primary_blue'])
    ax.plot(phi['median_value'].values / 0.4558, w['median_value'], 'o',
            markerfacecolor='w', markeredgecolor=cor['primary_blue'], markeredgewidth=1)

c_mapper = {'95%': cor['pale_black'],
            '75%': cor['light_black'],
            '25%': cor['primary_black'],
            'median': cor['black']}
for g, d in theo[theo['quantity'] == 'width_theory'].groupby('interval', sort=False):
    ax.fill_between(d['phiRb']/0.4558, d['lower'],
                    d['upper'], color=c_mapper[g], alpha=0.25)

ax.set_xlim([0.05/0.4558, 0.65])
ax.set_xlabel('RNA / protein', fontsize=6)
ax.set_ylabel('width [µm]', fontsize=6)
plt.savefig('./updated_width_theory.pdf')

# %%

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
mem = ms_data[ms_data['localization'] == 'membrane']
mem['rho'] = ms_data['mass_fg'].values / (2 * ms_data['surface_area'])
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_ylim([0, 5])
for g, d in mem.groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['rho'], **fmt)

for g, d in wt_data.groupby('carbon_source'):
    lam = d[d['quantity'] == 'growth_rate']
    rho_mem = d[d['quantity'] == 'rho_mem']

    ax.vlines(lam['median_value'].values, rho_mem['2.5%'],
              rho_mem['97.5%'], linewidth=1, color=cor['primary_blue'])
    ax.hlines(rho_mem['median_value'].values, lam['2.5%'],
              lam['97.5%'], linewidth=1, color=cor['primary_blue'])
    ax.plot(lam['median_value'].values, rho_mem['median_value'], 'o',
            markerfacecolor='w', markeredgecolor=cor['primary_blue'], markeredgewidth=1)

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('membrane protein density [fg/µm$^2$]', fontsize=6)

plt.savefig('./updated_membrane_plot.pdf')
