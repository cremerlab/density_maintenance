# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
ppc_cmap = {'95%': cor['pale_blue'], '75%': cor['light_blue'], '25%': cor['primary_blue'],
            'median': cor['blue']}
err_widths = {'95%': 0.5, '75%': 1, '25%': 1.5}

# Load the various dataset
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
shape_kdes = pd.read_csv(
    '../../data/mcmc/perturbation_shape_posterior_kde.csv')
param_kdes = pd.read_csv(
    '../../data/mcmc/perturbation_parameter_posterior_kde.csv')
percs = pd.read_csv('../../data/mcmc/perturbation_parameter_percentiles.csv')
_percs = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
percs = pd.concat([percs, _percs])
model = pd.read_csv('../../data/mcmc/literature_model_params.csv')
model = model[model['model'] == 'const_phi_mem']

# %%
# Set a utility function for plotting literature data and wildtypes


def canvas():
    fig, ax = plt.subplots(1, 4, figsize=(6, 1.5), sharex=True)
    axes = {'w_rep': ax[0], 'alpha_rep': ax[1],
            'phi_peri_rep': ax[2], 'm_peri_rep': ax[3]}
    for g, d in model[model['quantity'].isin(['w_rep', 'alpha_rep',
                                              'phi_peri_rep', 'm_peri_rep']) &
                      model['interval'].isin(ppc_cmap.keys())
                      ].groupby(['interval', 'quantity'], sort=False):
        axes[g[1]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                color=ppc_cmap[g[0]], alpha=0.25, zorder=1)
        if g[0] == 'median':
            axes[g[1]].plot(d['growth_rate_hr'], d['lower'], lw=1, color=ppc_cmap[g[0]],
                            alpha=0.25)

    for g, d in size_data.groupby(['source']):
        for i, p in enumerate(['width_um', 'aspect_ratio']):
            ax[i].plot(d['growth_rate_hr'], d[p], mapper[g]['m'], ms=3,
                       markeredgecolor=cor['primary_black'],
                       markerfacecolor=mapper[g]['c'], alpha=0.35, zorder=100)

    for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby(['dataset_name', 'condition', 'growth_rate_hr']):
        for i, p in enumerate(['mass_frac', 'mass_fg']):
            ax[i+2].plot(d['growth_rate_hr'], d[p], mapper[g[0]]['m'], ms=3,
                         markeredgecolor=cor['primary_black'],
                         markerfacecolor=mapper[g[0]]['c'], alpha=0.3, zorder=100)

    # Plot the wildtype inference percentiles
    for g, d in percs[(percs['strain'] == 'wildtype') & (percs['overexpression'] == 'none') &
                      (percs['inducer_conc'] == 0)].groupby(['carbon_source']):
        for i, p in enumerate(['width_mu', 'alpha_mu', 'phi_peri', 'm_peri']):
            _d = d[d['quantity'] == p]
            _g = d[d['quantity'] == 'growth_mu']
            med_lam = _g[_g['interval'] == 'median']['lower']
            med_p = _d[_d['interval'] == 'median']['lower']
            if (len(med_lam) == 0) | (len(med_p) == 0):
                continue
            for __g, __d in _d[_d['interval'].isin(err_widths.keys())].groupby(['interval'], sort=False):
                __lam = _g[_g['interval'] == __g]
                ax[i].vlines(med_lam, __d['lower'], __d['upper'], lw=err_widths[__g],
                             color=cor['blue'], zorder=998)
                ax[i].hlines(med_p, __lam['lower'], __lam['upper'], lw=err_widths[__g],
                             color=cor['blue'], zorder=998)
            ax[i].plot(med_lam, med_p, 'o', markerfacecolor='w', markeredgecolor=cor['blue'],
                       markeredgewidth=0.5, ms=3, zorder=999)

    for a in ax:
        a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
        a.set_xlim([0.15, 1.25])
    ax[0].set_title('average cell width', fontsize=6)
    ax[0].set_ylabel('$w$ [µm]', fontsize=6)
    ax[1].set_title('average aspect ratio', fontsize=6)
    ax[1].set_ylabel(r'$\alpha$', fontsize=6)
    ax[2].set_title('periplasmic allocation', fontsize=6)
    ax[2].set_ylabel(r'$\phi_{peri}$', fontsize=6)
    ax[3].set_title('periplasmic protein mass', fontsize=6)
    ax[3].set_ylabel(r'$m_{peri}$', fontsize=6)
    ax[0].set_ylim([0.45, 0.9])
    ax[1].set_ylim([1, 6])
    ax[2].set_ylim([0, 0.12])
    ax[3].set_ylim([0, 25])
    return fig, ax


# %%
delta3 = percs[(percs['strain'] == 'malE-rbsB-fliC-KO')
               & (percs['overexpression'] == 'none')]

fig, ax = canvas()
for g, d in delta3.groupby(['carbon_source']):
    lam = d[d['quantity'] == 'growth_mu']
    med_lam = lam[lam['interval'] == 'median']['lower']
    for i, p in enumerate(['width_mu', 'alpha_mu', 'phi_peri', 'm_peri']):
        med_p = d[(d['quantity'] == p) & (d['interval'] == 'median')]['lower']
        for _g, _d in d[d['interval'].isin(['95%', '75%', '25%']) &
                        (d['quantity'] == p)].groupby(['interval'], sort=False):
            _lam = lam[lam['interval'] == _g]
            ax[i].vlines(med_lam, _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=cor['primary_purple'], zorder=1000)
            ax[i].hlines(med_p, _lam['lower'], _lam['upper'], lw=err_widths[_g],
                         color=cor['primary_purple'], zorder=1000)
        ax[i].plot(med_lam, med_p, 'o', ms=3.5, color='w', markeredgecolor=cor['primary_purple'],
                   markeredgewidth=0.5, zorder=1001)
plt.tight_layout()
plt.savefig('../../figures/fig5_KO_comparison.pdf')

delta3 = percs[(percs['strain'] == 'malE-rbsB-fliC-KO')
               & (percs['overexpression'] == 'none')]

# %%
# MalE overexpression
overex = percs[(percs['strain'] == 'malE-rbsB-fliC-KO')]

fig, ax = canvas()
for g, d in overex.groupby(['carbon_source', 'inducer_conc', 'overexpression']):
    lam = d[d['quantity'] == 'growth_mu']
    med_lam = lam[lam['interval'] == 'median']['lower']
    if g[-1] == 'malE':
        c = cor['primary_green']
    else:
        c = cor['dark_green']

    for i, p in enumerate(['width_mu', 'alpha_mu', 'phi_peri', 'm_peri']):
        med_p = d[(d['quantity'] == p) & (d['interval'] == 'median')]['lower']
        for _g, _d in d[d['interval'].isin(['95%', '75%', '25%']) &
                        (d['quantity'] == p)].groupby(['interval'], sort=False):
            _lam = lam[lam['interval'] == _g]
            ax[i].vlines(med_lam, _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=c, zorder=1000)
            ax[i].hlines(med_p, _lam['lower'], _lam['upper'], lw=err_widths[_g],
                         color=c, zorder=1000)
        ax[i].plot(med_lam, med_p, 'o', ms=3.5, color='w', markeredgecolor=c,
                   markeredgewidth=0.5, zorder=1001)
plt.tight_layout()
plt.savefig('../../figures/fig5_OE_comparison.pdf')

# %%
# Overexpression
overex = percs[(percs['strain'] == 'wildtype') &
               (percs['overexpression'] == 'lacZ')]
fig, ax = canvas()
for g, d in overex.groupby(['carbon_source', 'inducer_conc', 'overexpression']):
    lam = d[d['quantity'] == 'growth_mu']
    med_lam = lam[lam['interval'] == 'median']['lower']
    c = cor['primary_black']
    for i, p in enumerate(['width_mu', 'alpha_mu', 'phi_peri', 'm_peri']):
        med_p = d[(d['quantity'] == p) & (d['interval'] == 'median')]['lower']
        for _g, _d in d[d['interval'].isin(['95%', '75%', '25%']) &
                        (d['quantity'] == p)].groupby(['interval'], sort=False):
            _lam = lam[lam['interval'] == _g]
            ax[i].vlines(med_lam, _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=c, zorder=1000)
            ax[i].hlines(med_p, _lam['lower'], _lam['upper'], lw=err_widths[_g],
                         color=c, zorder=1000)
        ax[i].plot(med_lam, med_p, 'o', ms=3.5, color='w', markeredgecolor=c,
                   markeredgewidth=0.5, zorder=1001)
plt.tight_layout()
plt.savefig('../../figures/fig5_lacZ_comparison.pdf')

# %%
# # Overexpression
# lpp = percs[(percs['strain'] == 'lpp14')]
# fig, ax = canvas()
# for g, d in lpp.groupby(['carbon_source']):
#     if g == 'LB':
#         continue
#     lam = d[d['quantity'] == 'growth_mu']

#     med_lam = lam[lam['interval'] == 'median']['lower']
#     c = cor['primary_red']
#     for i, p in enumerate(['width_mu', 'alpha_mu', 'phi_peri', 'm_peri']):
#         med_p = d[(d['quantity'] == p) & (d['interval'] == 'median')]['lower']
#         for _g, _d in d[d['interval'].isin(['95%', '75%', '25%']) &
#                         (d['quantity'] == p)].groupby(['interval'], sort=False):
#             _lam = lam[lam['interval'] == _g]
#             ax[i].vlines(med_lam, _d['lower'], _d['upper'], lw=err_widths[_g],
#                          color=c, zorder=1000)
#             ax[i].hlines(med_p, _lam['lower'], _lam['upper'], lw=err_widths[_g],
#                          color=c, zorder=1000)
#         ax[i].plot(med_lam, med_p, 'o', ms=3.5, color='w', markeredgecolor=c,
#                    markeredgewidth=0.5, zorder=1001)
# plt.tight_layout()
# plt.savefig('../../figures/fig5_lpp_comparison.pdf')
