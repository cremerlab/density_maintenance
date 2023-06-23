# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# %%
# Add our data with 95%
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[~size_data['source'].isin(
    ['Si et al. 2017', 'Taher-Araghi et al. 2015'])]
wt_data = pd.read_csv('../../data/mcmc/perturbation_parameter_percentiles.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
_wt_data = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
wt_data = pd.concat([wt_data, _wt_data])
_data = pd.concat([wt_data, _wt_data])
wt_data = wt_data[(wt_data['strain'] == 'wildtype') &
                  (wt_data['overexpression'] == 'none') &
                  (wt_data['inducer_conc'] == 0) &
                  (wt_data['interval'].isin(['median', '95%']))]
model = pd.read_csv('../../data/mcmc/literature_model_params.csv')
model = model[(model['volume_scale'] == 'linear_width')
              & (model['model'] == 'const_phi_mem')]

# %%
width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
size_data['mass_fraction'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width_um'] = width_popt[1] + \
    width_popt[0] * phiRb_data['growth_rate_hr']
ms = ms_data[ms_data['localization'] == 'ribosomal sector'].copy()
ms.rename(columns={'dataset_name': 'source', 'width': 'width_um',
          'mass_frac': 'mass_fraction'}, inplace=True)
merged = pd.concat([size_data[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']],
                   phiRb_data[['source', 'width_um',
                               'mass_fraction', 'growth_rate_hr']],
                   ms[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']]])
merged = merged[(merged['source'] != 'Si et al. 2017') & (
    merged['source'] != 'Taheri-Araghi et al. 2015')]
size_data = size_data[size_data['source'] != 'Si et al. 2017']
size_data = size_data[size_data['source'] != 'Taheri-Araghi et al. 2015']
# %%
popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], np.log(size_data['volume_um3']))
lam_range = np.linspace(0, 2.5, 200)
fit = np.exp(popt[1] + popt[0] * lam_range)
fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.75))

for g, d in size_data.groupby('source'):
    ax.plot(d['growth_rate_hr'], d['volume_um3'], mapper[g]['m'], ms=5,
            markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
            color=mapper[g]['c'], alpha=0.75, label=g)

for g, d in wt_data[wt_data['quantity'].isin(['volume_mu', 'growth_mu'])].groupby(['carbon_source']):
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]
    vol_med = med[med['quantity'] == 'volume_mu']['lower'].values[0]
    vol_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'volume_mu')]
    ax.hlines(vol_med, lam_95['lower'], lam_95['upper'], lw=1,
              color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax.vlines(lam_med, vol_95['lower'], vol_95['upper'], lw=1,
              color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax.plot(lam_med, vol_med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
            ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)

ax.errorbar([], [], xerr=[], yerr=[], fmt='o', color=cor['primary_blue'], markeredgecolor=cor['primary_blue'],
            markeredgewidth=1, markerfacecolor='w', ms=4.5, capsize=0, label='our data')
ax.plot(lam_range, fit, lw=1.5, color=cor['primary_black'], zorder=10)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('average cell volume [µm$^{3}$]', fontsize=6)
ax.legend(fontsize=4)  # bbox_to_anchor=(1, 1), handletextpad=0.1)
ax.set_ylim([-0.1, 6])
ax.set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
plt.savefig('../../figures/figXG_volume_plot_our_data.pdf',
            bbox_inches='tight')
plt.savefig('../../figures/figXG_volume_plot_our_data.png',
            bbox_inches='tight', dpi=300)

# %%
fig, ax = plt.subplots(1, 2, figsize=(3.25, 1.5), sharex=True)
ax[0].set_ylim([0, 0.15])
ax[0].set_xlim([0, 2])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2])
ax[1].set_ylim([0, 0.15])
ax[1].set_ylabel('proteome mass fraction', fontsize=6)

for g, d in ms_data[ms_data['localization'].isin(['inner membrane', 'outer membrane', 'periplasm']
                                                 )].groupby(['dataset_name', 'localization']):
    print(g)
    if g[1] == 'inner membrane':
        c = cor['pale_blue']
        a = ax[0]
    elif g[1] == 'outer membrane':
        c = cor['dark_blue']
        a = ax[0]
    else:
        c = cor['purple']
        a = ax[1]
    a.plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], markeredgecolor=cor['primary_black'],
           markerfacecolor=c, alpha=0.5, markeredgewidth=0.5, ms=4.5, label='__nolegend__')

for g, d in ms_data.groupby(['dataset_name']):
    ax[1].plot([], [], mapper[g]['m'], markeredgewidth=0.5, markerfacecolor=mapper[g]['c'],
               alpha=0.5, ms=4.5, label=g, markeredgecolor=cor['primary_black'])


for g, d in wt_data[wt_data['quantity'].isin(['phi_peri', 'growth_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]
    phi_med = med[med['quantity'] == 'phi_peri']['lower'].values[0]
    phi_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'phi_peri')]
    ax[1].hlines(phi_med, lam_95['lower'], lam_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].vlines(lam_med, phi_95['lower'], phi_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].plot(lam_med, phi_med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
               ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)

ax[1].errorbar([], [], xerr=[], yerr=[], fmt='o', color=cor['primary_blue'], markeredgecolor=cor['primary_blue'],
               markeredgewidth=1, markerfacecolor='w', ms=4.5, capsize=0, label='our data')
ax[0].plot([], [], 's', color=cor['pale_blue'], label='inner membrane', ms=4)
ax[0].plot([], [], 's', color=cor['dark_blue'], label='outer membrane', ms=4)
ax[0].legend(fontsize=4, handletextpad=0.1)
ax[1].legend(fontsize=4, handletextpad=0.1)


for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('proteome mass fraction', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/FigXG_allocation_nomodel.pdf', bbox_inches='tight')

# %%
# %%
fig, ax = plt.subplots(1, 2, figsize=(3.25, 1.5), sharex=True)
ax[0].set_ylim([0, 0.25])
ax[0].set_xlim([0, 2])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2])
ax[1].set_ylim([0, 0.15])
ax[1].set_ylabel('proteome mass fraction', fontsize=6)


inter_colors = {'95%': cor['pale_black'],
                '75%': cor['light_black'], '25%': cor['primary_black'],
                'median': cor['black']}

for g, d in model[model['quantity'].isin(['phi_mem_rep', 'phi_peri_rep']) &
                  model['interval'].isin(inter_colors.keys())
                  ].groupby(['quantity', 'interval'], sort=False):

    if g[0] == 'phi_mem_rep':
        a = ax[0]
    else:
        a = ax[1]

    if g[1] != 'median':
        a.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=inter_colors[g[1]],
                       alpha=0.35)
    else:
        a.plot(d['growth_rate_hr'], d['lower'], '-',
               color=inter_colors[g[1]], alpha=0.35)

for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm']
                                                 )].groupby(['dataset_name', 'localization']):
    c = mapper[g[0]]['c']
    if g[1] == 'membrane':
        a = ax[0]
    else:
        a = ax[1]
    a.plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], markeredgecolor=cor['primary_black'],
           markerfacecolor=c, alpha=0.5, markeredgewidth=0.5, ms=4.5, label='__nolegend__')
for g, d in ms_data.groupby(['dataset_name']):
    ax[1].plot([], [], mapper[g]['m'], markeredgewidth=0.5, markerfacecolor=mapper[g]['c'],
               alpha=0.5, ms=4.5, label=g, markeredgecolor=cor['primary_black'])


for g, d in wt_data[wt_data['quantity'].isin(['phi_peri', 'growth_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]
    phi_med = med[med['quantity'] == 'phi_peri']['lower'].values[0]
    phi_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'phi_peri')]
    ax[1].hlines(phi_med, lam_95['lower'], lam_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].vlines(lam_med, phi_95['lower'], phi_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].plot(lam_med, phi_med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
               ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)

ax[1].errorbar([], [], xerr=[], yerr=[], fmt='o', color=cor['primary_blue'], markeredgecolor=cor['primary_blue'],
               markeredgewidth=1, markerfacecolor='w', ms=4.5, capsize=0, label='our data')
ax[1].legend(fontsize=4, handletextpad=0.1)


for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('proteome mass fraction', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/FigXG_allocation_model.pdf', bbox_inches='tight')


# %%
fig, ax = plt.subplots(2, 2, figsize=(3, 1.75), sharex=True)
ax[0, 0].set_xlim([-0.1, 2])
ax[0, 0].set_xticks([0, 0.5, 1,  1.5, 2, 2])
ax[0, 0].set_ylim([0.3, 1.3])
ax[0, 0].set_yticks([0.5, 0.75, 1.0])
ax[1, 0].set_ylim([0.5, 6])
ax[1, 0].set_yticks([2, 4, 6])
ax[1, 1].set_ylim([1, 8])
ax[0, 1].set_ylim([0, 11])
ax[1, 1].set_yticks([1, 3, 5, 7])
ax[0, 1].set_yticks([1, 5, 9])
ax[1, 0].set_yticklabels(['2.00', '4.00', '6.00'])
ax[1, 0].set_ylabel(r'$\langle \ell \rangle$ [µm]', fontsize=6)
ax[0, 0].set_ylabel(r'$\langle w \rangle$ [µm]', fontsize=6)
ax[1, 0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1, 1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0, 1].set_ylabel('S/V [µm$^{-1}$]', fontsize=6)
ax[1, 1].set_ylabel('avg. aspect ratio', fontsize=6)

axes = {'ell_rep': ax[1, 0], 'w_rep': ax[0, 0],
        'sav_rep': ax[0, 1], 'alpha_rep': ax[1, 1]}
for g, d in model[model['quantity'].isin(['ell_rep', 'w_rep', 'sav_rep', 'alpha_rep']) &
                  model['interval'].isin(inter_colors.keys())
                  ].groupby(['quantity', 'interval'], sort=False):

    a = axes[g[0]]
    if g[1] != 'median':
        a.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=inter_colors[g[1]],
                       alpha=0.35)
    else:
        a.plot(d['growth_rate_hr'], d['lower'], '-',
               color=inter_colors[g[1]], alpha=0.35)

for g, d in size_data.groupby(['source']):
    ax[0, 0].plot(d['growth_rate_hr'], d['width_um'], mapper[g]['m'],
                  markeredgecolor=cor['primary_black'],
                  markerfacecolor=mapper[g]['c'], ms=4.5, alpha=0.5)
    ax[1, 0].plot(d['growth_rate_hr'], d['length_um'], mapper[g]['m'],
                  markeredgecolor=cor['primary_black'],
                  markerfacecolor=mapper[g]['c'], ms=4.5, alpha=0.5)
    ax[0, 1].plot(d['growth_rate_hr'], d['surface_to_volume'], mapper[g]['m'],
                  markeredgecolor=cor['primary_black'],
                  markerfacecolor=mapper[g]['c'], ms=4.5, alpha=0.5)
    ax[1, 1].plot(d['growth_rate_hr'], d['length_um'] / d['width_um'], mapper[g]['m'],
                  markeredgecolor=cor['primary_black'],
                  markerfacecolor=mapper[g]['c'], ms=4.5, alpha=0.5)

vars = {'width_mu': ax[0, 0], 'length_mu': ax[1, 0],
        'sav_mu': ax[0, 1], 'alpha_mu': ax[1, 1]}

for g, d in wt_data[wt_data['quantity'].isin(['width_mu', 'length_mu', 'sav_mu', 'alpha_mu', 'growth_mu'])].groupby(['carbon_source']):
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]

    for _v, a in vars.items():
        _med = med[med['quantity'] == _v]['lower'].values[0]
        _95 = d[(d['interval'] == '95%') & (d['quantity'] == _v)]

        a.hlines(_med, lam_95['lower'], lam_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
        a.vlines(lam_med, _95['lower'], _95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
        a.plot(lam_med, _med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
               ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)

plt.savefig('../../figures/FigXG_sizes.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3.25, 1.75), sharex=True)
ax[0].set_ylim([0, 0.15])
ax[0].set_xlim([0, 2])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2])

ax[1].set_ylim([0, 0.15])
ax[1].set_ylabel('proteome mass fraction', fontsize=6)

for g, d in ms_data[ms_data['localization'].isin(['inner membrane', 'outer membrane', 'periplasm']
                                                 )].groupby(['dataset_name', 'localization']):
    print(g)
    if g[1] == 'inner membrane':
        c = cor['pale_blue']
        a = ax[0]
    elif g[1] == 'outer membrane':
        c = cor['dark_blue']
        a = ax[0]
    else:
        c = cor['purple']
        a = ax[1]
    a.plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], markeredgecolor=cor['primary_black'],
           markerfacecolor=c, alpha=0.5, markeredgewidth=0.5, ms=4.5, label='__nolegend__')

for g, d in ms_data.groupby(['dataset_name']):
    ax[1].plot([], [], mapper[g]['m'], markeredgewidth=0.5, markerfacecolor=mapper[g]['c'],
               alpha=0.5, ms=4.5, label=g, markeredgecolor=cor['primary_black'])


for g, d in wt_data[wt_data['quantity'].isin(['m_peri', 'growth_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]
    phi_med = med[med['quantity'] == 'phi_peri']['lower'].values[0]
    phi_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'phi_peri')]
    ax[1].hlines(phi_med, lam_95['lower'], lam_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].vlines(lam_med, phi_95['lower'], phi_95['upper'], lw=1,
                 color=cor['primary_blue'], label='__nolegend__', zorder=1000)
    ax[1].plot(lam_med, phi_med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
               ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)

ax[1].errorbar([], [], xerr=[], yerr=[], fmt='o', color=cor['primary_blue'], markeredgecolor=cor['primary_blue'],
               markeredgewidth=1, markerfacecolor='w', ms=4.5, capsize=0, label='our data')
ax[0].plot([], [], 's', color=cor['pale_blue'], label='inner membrane', ms=4)
ax[0].plot([], [], 's', color=cor['dark_blue'], label='outer membrane', ms=4)
ax[0].legend(fontsize=4, handletextpad=0.1)
ax[1].legend(fontsize=4, handletextpad=0.1)


for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('proteome mass fraction', fontsize=6)
plt.tight_layout()
# plt.savefig('../../figures/FigXG_allocation_nomodel.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)
ax[0].set_xlim([0, 2.1])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('membrane protein density\n[fg / µm$^2$]', fontsize=6)
ax[1].set_ylabel('periplasmic protein density\n[fg / µm$^3$]', fontsize=6)


for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm'])].groupby(['dataset_name', 'localization']):
    if g[1] == 'membrane':
        v = d['mass_fg'] / (2 * d['surface_area'])
        a = ax[0]
    else:
        v = d['mass_fg'] / (0.0246 * d['surface_area'])
        a = ax[1]

    a.plot(d['growth_rate_hr'], v, mapper[g[0]]['m'],
           markeredgecolor=cor['primary_black'], alpha=0.5, markerfacecolor=mapper[g[0]]['c'],
           ms=4.5)


for g, d in model[model['quantity'].isin(['rho_mem_rep', 'rho_peri_rep']) &
                  model['interval'].isin(inter_colors.keys())
                  ].groupby(['quantity', 'interval'], sort=False):

    if g[0] == 'rho_mem_rep':
        a = ax[0]
    else:
        a = ax[1]

    if g[1] != 'median':
        a.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=inter_colors[g[1]],
                       alpha=0.35)
    else:
        a.plot(d['growth_rate_hr'], d['lower'], '-',
               color=inter_colors[g[1]], alpha=0.35)


for g, d in wt_data[wt_data['quantity'].isin(['rho_peri', 'rho_mem_rep', 'growth_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    med = d[d['interval'] == 'median']
    lam_med = med[med['quantity'] == 'growth_mu']['lower'].values[0]
    lam_95 = d[(d['interval'] == '95%') & (d['quantity'] == 'growth_mu')]

    for i, _v in enumerate(['rho_peri']):
        _med = med[med['quantity'] == _v]['lower'].values[0]
        _95 = d[(d['interval'] == '95%') & (d['quantity'] == _v)]
        ax[1].hlines(_med, lam_95['lower'], lam_95['upper'], lw=1,
                     color=cor['primary_blue'], label='__nolegend__', zorder=1000)
        ax[1].vlines(lam_med, _95['lower'], _95['upper'], lw=1,
                     color=cor['primary_blue'], label='__nolegend__', zorder=1000)
        ax[1].plot(lam_med, _med, 'o', markeredgecolor=cor['primary_blue'], markerfacecolor='w',
                   ms=4.5, label='__nolegend__', zorder=1000, markeredgewidth=1)


ax[0].set_ylim([0, 6])
ax[1].set_ylim([0, 150])
plt.tight_layout()
plt.savefig('../../figures/FigXG_densities_model.pdf', bbox_inches='tight')

# %%
# Bump chart
fig, ax = plt.subplots(2, 2, figsize=(2.5, 1.5), sharex=True)

wt_ace = wt_data[wt_data['carbon_source'] == 'acetate']
d3_ace = _data[(_data['strain'] == 'malE-rbsB-fliC-KO') &
               (_data['overexpression'] == 'none') &
               (_data['inducer_conc'] == 0) &
               (_data['carbon_source'] == 'acetate')]
malE_ace = _data[(_data['strain'] == 'malE-rbsB-fliC-KO') &
                 (_data['overexpression'] == 'malE') &
                 (_data['inducer_conc'] == 100) &
                 (_data['carbon_source'] == 'acetate')]
lpp_ace = _data[(_data['strain'] == 'lpp14') &
                (_data['overexpression'] == 'none') &
                (_data['inducer_conc'] == 0) &
                (_data['carbon_source'] == 'acetate')]

axes = {'width_mu': ax[0, 0], 'length_mu': ax[0, 1],
        'alpha_mu': ax[1, 0], 'm_peri': ax[1, 1]}
meds = {}
for i, (_d, c) in enumerate(zip([wt_ace, lpp_ace, d3_ace, malE_ace],
                                [cor['primary_blue'], cor['primary_red'], cor['primary_purple'],
                                 cor['primary_green']])):
    for g, d in _d[_d['quantity'].isin(axes.keys())].groupby(['quantity']):
        _med = d[d['interval'] == 'median']
        _95 = d[d['interval'] == '95%']
        if i == 0:
            meds[g] = _med['lower'].values[0]
        axes[g].vlines(i, (_95['lower'] - meds[g])/meds[g],
                       (_95['upper'] - meds[g])/meds[g], lw=0.5, color=c)
        axes[g].plot(i, (_med['lower'] - meds[g])/meds[g], 'o', markeredgecolor=c, markerfacecolor='w',
                     ms=4.5)

for a in ax.ravel()[:-1]:
    a.set_ylim([-0.15, 0.25])
    # a.set_ylabel("% change from WT", fontsize=6)
    a.set_yticks([-0.1, 0, 0.1, 0.2])
    a.set_yticklabels(['-10%', '0%', '+10%', '+20%'])
    a.set_xlim([-0.2, 3.2])
ax[1, 1].set_yticks([-0.5, -0.25, 0])
ax[1, 1].set_yticklabels(['-50%', '-25%', '0%'])
ax[0, 0].set_title('cell width', fontsize=6)
ax[0, 1].set_title('cell length', fontsize=6)
ax[1, 0].set_title('cell aspect ratio', fontsize=6)
ax[1, 1].set_title('periplasmic protein mass', fontsize=6)

ax[1, 0].set_xticks([0, 1, 2, 3])
ax[1, 1].set_xticks([0, 1, 2, 3])
ax[1, 0].set_xticklabels(
    ['WT', 'lpp$^{14}$', '$\Delta$3',  '$\Delta$3\n+ malE'])
ax[1, 1].set_xticklabels(
    ['WT', 'lpp$^{14}$', '$\Delta$3',  '$\Delta$3\n+ malE'])
plt.subplots_adjust(hspace=0.2, wspace=0.35)
plt.savefig('../../figures/FigXG_perc_change_sizes.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_ylim([0, 25])
ax.set_xlim([0, 2])
for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['mass_fg'], mapper[g]['m'],
            markeredgecolor=cor['primary_black'], markerfacecolor=mapper[g]['c'],
            alpha=0.5, ms=4.5, label='__nolegend__')

c = {'wildtype, none': cor['primary_blue'],
     'malE-rbsB-fliC-KO, none': cor['primary_purple'],
     'malE-rbsB-fliC-KO, malE': cor['primary_green'],
     'malE-rbsB-fliC-KO, rbsB': cor['primary_gold'],
     'lpp14, none': cor['primary_red']}

for g, d in model[(model['quantity'] == 'm_peri_rep')
                  & model['interval'].isin(inter_colors.keys())].groupby(['interval'], sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                    color=inter_colors[g], alpha=0.5, label='__nolegend__')

labels = []
for g, d in _data[_data['quantity'].isin(['m_peri', 'growth_mu']) &
                  _data['interval'].isin(['median', '95%']) &
                  _data['overexpression'].isin(['none', 'malE', 'rbsB'])
                  ].groupby(['strain', 'carbon_source',
                             'overexpression', 'inducer_conc']):
    if g[1] == 'LB':
        continue
    _c = c[f'{g[0]}, {g[2]}']

    if g[0] == 'malE-rbsB-fliC-KO':
        label = f'$\Delta$3'
    else:
        label = f'{g[0]}'
    if g[2] != 'none':
        label += f' + {g[2]}'
    if label not in labels:
        labels.append(label)
    else:
        label = '__nolegend__'
    _med = d[d['interval'] == 'median']
    _95 = d[d['interval'] == '95%']
    ax.vlines(_med[_med['quantity'] == 'growth_mu']['lower'].values[0],
              _95[_95['quantity'] == 'm_peri']['lower'].values[0],
              _95[_95['quantity'] == 'm_peri']['upper'].values[0],
              color=_c, lw=0.75, label='__nolegend__')
    ax.hlines(_med[_med['quantity'] == 'm_peri']['lower'].values[0],
              _95[_95['quantity'] == 'growth_mu']['lower'].values[0],
              _95[_95['quantity'] == 'growth_mu']['upper'].values[0],
              color=_c, lw=0.75, label='__nolegend__')
    ax.plot(_med[_med['quantity'] == 'growth_mu']['lower'].values[0],
            _med[_med['quantity'] == 'm_peri']['lower'].values[0],
            'o', ms=4.5, markeredgecolor=_c, markerfacecolor='w',
            markeredgewidth=0.75, label=label)

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('periplasmic protein [fg / cell]', fontsize=6)
ax.legend(fontsize=4, loc='upper right')
plt.savefig('../../figures/FigXG_m_peri.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 1, figsize=(1.25, 1.9), sharex=True, sharey=True)
ax.set_xticks([0, 1])
ax.set_ylim([-10, 13])
ax.set_yticks([-10, -5, 0, 5, 10])
ax.set_ylabel('log$_2$ enrichment in supernatant', fontsize=6)
ax.set_xticklabels(['other\nproteins', 'periplasmic\nproteins'], fontsize=6)

colors = ['black', 'purple']
ms_conf['sublocalization'] = 'periplasmic'
ms_conf.loc[ms_conf['localization'] != 'periplasm',
            'sublocalization'] = 'non-periplasmic'
for i, (g, d) in enumerate(ms_conf.groupby(['sublocalization'])):
    print(g)
    _ = ax.plot(np.random.normal(0, 0.1, len(d)) + i, np.log2(d['relative_signal']), '.', color=cor[f'primary_{colors[i]}'], markeredgewidth=0, alpha=0.75,
                ms=2.5)

    # for i, p in enumerate(['peri_enrichment', 'cyto_enrichment']):
    # _ = axes[g].plot(np.log2(d[p]), np.random.normal(
    # 0, 0.1, len(d)) + i, '.', color=cor[f'primary_{colors[i]}'], markeredgewidth=0, alpha=0.5,
    # ms=3)
ax.set_title('mass spectrometry\nconfirmation', fontsize=6)
# plt.tight_layout()
plt.savefig('../../figures/FigSX_MS_periplasm_validation_plots.pdf')

# %%
