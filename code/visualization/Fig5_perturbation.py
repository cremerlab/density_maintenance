# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
err_widths = {'95%': 0.5, '75%': 1, '25%': 1.5}
pert_cors = {'wildtype': {'lacZ': cor['dark_black'],
                          'none': cor['primary_blue']},
             'malE-rbsB-fliC-KO': {'none': cor['primary_purple'],
                                   'malE': cor['primary_green'],
                                   'rbsB': cor['primary_gold']},
             'lpp14': {'none': cor['primary_red']}}
perc_cors = {'95%': cor['pale_black'], '75%': cor['light_black'],
             '25%': cor['primary_black'], 'median': cor['black']}
# Load datasets
growth_params = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
model_params = pd.read_csv(
    '../../data/mcmc/perturbation_parameter_percentiles.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
kdes = pd.read_csv('../../data/mcmc/perturbation_shape_posterior_kde.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
ppcs = ppcs[ppcs['model'] == 'const_phi_mem']
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']


# %%
# Size KDEs
fig, ax = plt.subplots(3, 1, figsize=(2, 3), sharex=True)
axes = {'width_mu': ax[0], 'length_mu': ax[1], 'alpha_mu': ax[2]}
for a in ax:
    a.set_xlim([0.6, 1.3])
    a.set_ylim([-0.1, 1.25])
wt = kdes[(kdes['strain'] == 'wildtype') & (kdes['carbon_source'] == 'acetate') &
          (kdes['overexpression'] == 'none')]

deltas = []
for i, (g, d) in enumerate(kdes[(kdes['carbon_source'] == 'acetate') &
                                kdes['inducer_conc'].isin([0, 50])
                                ].groupby(['strain', 'overexpression'])):
    if g[1] != 'none':
        d = d[d['inducer_conc'] > 0]
    for _g, _d in d[d['parameter'].isin(list(axes.keys()))].groupby(['parameter']):
        _wt = wt[wt['parameter'] == _g]
        _wt = _wt.iloc[np.argmax(_wt['kde'])]['value']
        if (g[0] == 'lpp14'):
            maxdiff = 1 - (_d.iloc[np.argmax(_d['kde'].values)]['value'] / _wt)
            print(_g, maxdiff)
            deltas.append(maxdiff)

        elif (g[0] == 'malE-rbsB-fliC-KO') & (g[1] == 'rbsB') & (_g == 'length_mu'):
            maxdiff = 1 - (_d.iloc[np.argmax(_d['kde'].values)]['value'] / _wt)
            print(_g, maxdiff)
            deltas.append(maxdiff)

        axes[_g].fill_between(_d.value/_wt, np.zeros(len(_d)), _d.kde/_d.kde.max(),
                              color=pert_cors[g[0]][g[1]], lw=1, alpha=0.25,
                              zorder=1001 + i)
        axes[_g].plot(_d.value/_wt, _d.kde/_d.kde.max(), '-',
                      color=pert_cors[g[0]][g[1]], lw=1, zorder=1000 + i)

for a in ax:
    a.set_yticks([])
    a.set_ylabel('$\propto$ probability', fontsize=6)


ax[2].set_xlabel(
    'fold-change from most probable\nwildtype value', fontsize=6)
ax[0].set_title('average cell width in acetate', fontsize=6)
ax[1].set_title('average cell length in acetate', fontsize=6)
ax[2].set_title('average cell aspect ratio in acetate', fontsize=6)
plt.savefig('../../figures/Fig5_perturbation_size_kdes.pdf')

# %%
fig, ax = plt.subplots(3, 1, figsize=(2, 3), sharex=True)
ax = ax.ravel()
for a in ax:
    a.set_xlim([0.15, 1.25])
ax[0].set_ylim([0.3, 1.3])
ax[1].set_ylim([1, 5])

ax[2].set_ylim([0.5, 8])
ax[0].set_ylabel('average\nwidth [µm]', fontsize=6)
ax[1].set_ylabel('average\nlength [µm]\n', fontsize=6)
ax[2].set_ylabel('average\naspect ratio\n', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

pert_cors = {'wildtype': {'lacZ': cor['dark_black'],
                          'none': cor['blue']},
             'malE-rbsB-fliC-KO': {'none': cor['primary_purple'],
                                   'malE': cor['primary_green'],
                                   'rbsB': cor['primary_gold']},
             'lpp14': {'none': cor['primary_red']}}

# Plot the percentiles
axes = {'w_rep': ax[0], 'ell_rep': ax[1], 'alpha_rep': ax[2]}
for i, (g, d) in enumerate(ppcs[ppcs['quantity'].isin(['w_rep', 'ell_rep', 'alpha_rep']) &
                                ppcs['interval'].isin(perc_cors.keys())
                                ].groupby(['quantity', 'interval'], sort=False)):
    if g[1] != 'median':
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], alpha=0.25,
                                color=perc_cors[g[1]], zorder=100)
    else:
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'], lw=1, alpha=0.25,
                        color=perc_cors[g[1]])

for g, d in size_data.groupby(['source']):
    for i, p in enumerate(['width_um', 'length_um', 'aspect_ratio']):
        ax[i].plot(d['growth_rate_hr'], d[p], mapper[g]['m'], ms=4,
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   color=mapper[g]['c'], alpha=0.45, zorder=500)

for g, d in growth_params.groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source', 'inducer', 'temperature'], sort=False):
    sizes = model_params[(model_params['strain'] == g[0]) &
                         (model_params['overexpression'] == g[1]) &
                         (model_params['inducer_conc'] == g[2]) &
                         (model_params['carbon_source'] == g[3]) &
                         (model_params['quantity'].isin(['width_mu', 'length_mu', 'alpha_mu']))]
    if g[1] == 'lacZ':
        continue
    med_growth = d[d['interval'] == 'median']
    for i, p in enumerate(['width_mu', 'length_mu', 'alpha_mu']):
        med_p = sizes[(sizes['quantity'] == p) &
                      (sizes['interval'] == 'median')]
        ax[i].plot(med_growth['lower'], med_p['lower'], 'o', ms=3, markeredgecolor=pert_cors[g[0]][g[1]],
                   markeredgewidth=0.5, markerfacecolor='white', zorder=1000)
        for _g, _d in sizes[(sizes['interval'] != 'median') & (sizes['quantity'] == p)].groupby(['interval']):
            ax[i].vlines(med_growth['lower'], _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=pert_cors[g[0]][g[1]], zorder=999)
        for _g, _d in d[d['interval'] != 'median'].groupby(['interval']):
            ax[i].hlines(med_p['lower'], _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=pert_cors[g[0]][g[1]], zorder=999)
plt.tight_layout()
plt.subplots_adjust(hspace=0.2)
plt.savefig('../../figures/Fig5_perturbation_dimensions.pdf',
            bbox_inches='tight')
# %%
# Protein Measurements
fig, ax = plt.subplots(1, 3, figsize=(7, 2), sharex=True)
ax[0].set_ylim([0, 25])
ax[1].set_ylim([0, 0.15])
ax[0].set_xlim([0.2, 1.3])
ax[1].set_xlim([0.2, 1.3])
# ax[2].set_xlim([0.5, 1.25])
ax[2].set_ylim([25, 175])
ax[0].set_ylabel('$m_{peri}$ [fg / cell]', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$', fontsize=6)
ax[2].set_ylabel(r'$\rho_{peri}$', fontsize=6)
ax[0].set_title('periplasmic protein mass', fontsize=6)
ax[1].set_title('allocation to periplasmic protein', fontsize=6)
ax[2].set_title('periplasmic protein density', fontsize=6)
for a in ax[:-1]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[-1].set_xlabel('average cell width [µm]', fontsize=6)

for i, (g, d) in enumerate(ppcs[ppcs['quantity'].isin(['m_peri_rep', 'phi_peri_rep', 'rho_peri_rep']) &
                                ppcs['interval'].isin(perc_cors.keys())
                                ].groupby(['quantity', 'interval'], sort=False)):
    x = d['growth_rate_hr']
    if g[0] == 'm_peri_rep':
        a = ax[0]
    elif g[0] == 'phi_peri_rep':
        a = ax[1]
    else:
        a = ax[2]
    if g[1] != 'median':
        a.fill_between(x, d['lower'], d['upper'], alpha=0.25,
                       color=perc_cors[g[1]], zorder=100)
    else:
        a.plot(x, d['lower'], lw=1, alpha=0.25,
               color=perc_cors[g[1]])


for g, d in ms_data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    peri = d[d['localization'] == 'periplasm']
    tot = d[d['localization'].isin(['cytoplasm', 'envelope'])]
    ax[0].plot(peri['growth_rate_hr'], peri['mass_fg'], mapper[g[0]]['m'],
               ms=4, color=mapper[g[0]]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.25, alpha=0.45)
    ax[1].plot(peri['growth_rate_hr'], peri['mass_fg']/tot['mass_fg'].sum(), mapper[g[0]]['m'],
               ms=4, color=mapper[g[0]]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.25, alpha=0.45)
    ax[2].plot(peri['growth_rate_hr'], peri['mass_fg']/(peri['surface_area'] * 0.0246), mapper[g[0]]['m'],
               ms=4, color=mapper[g[0]]['c'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.25, alpha=0.45)

for g, d in growth_params.groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source'], sort=False):
    pars = model_params[(model_params['strain'] == g[0]) &
                        (model_params['overexpression'] == g[1]) &
                        (model_params['inducer_conc'] == g[2]) &
                        (model_params['carbon_source'] == g[3]) &
                        (model_params['quantity'].isin(['m_peri', 'phi_peri', 'rho_peri']))]
    if len(pars) == 0:
        continue
    if g[1] == 'lacZ':
        continue
    med_growth = d[d['interval'] == 'median']
    for i, p in enumerate(['m_peri', 'phi_peri', 'rho_peri']):
        if p == 'm_peri_rep':
            prefactor = 1E9
        else:
            prefactor = 1
        med_p = pars[(pars['quantity'] == p) &
                     (pars['interval'] == 'median')]

        ax[i].plot(med_growth['lower'], prefactor * med_p['lower'], 'o', ms=3.5, markeredgecolor=pert_cors[g[0]][g[1]],
                   markeredgewidth=0.5, markerfacecolor='white', zorder=1000)
        for _g, _d in pars[(pars['interval'] != 'median') & (pars['quantity'] == p)].groupby(['interval']):
            ax[i].vlines(med_growth['lower'], prefactor * _d['lower'], prefactor * _d['upper'], lw=err_widths[_g],
                         color=pert_cors[g[0]][g[1]], zorder=999)
        for _g, _d in d[d['interval'] != 'median'].groupby(['interval']):
            ax[i].hlines(prefactor * med_p['lower'],  _d['lower'], _d['upper'], lw=err_widths[_g],
                         color=pert_cors[g[0]][g[1]], zorder=999)

# for g, d in model_params[model_params['quantity'] == 'width_mu'].groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc']):
#     pars = model_params[(model_params['strain'] == g[0]) &
#                         (model_params['overexpression'] == g[2]) &
#                         (model_params['inducer_conc'] == g[3]) &
#                         (model_params['carbon_source'] == g[1]) &
#                         (model_params['quantity'] == 'rho_peri')]
#     if len(pars) == 0:
#         continue
#     med_width = d[d['interval'] == 'median']
#     ax[2].plot(med_width['lower'], pars[(pars['interval'] == 'median')]['lower'], 'o', ms=3.5, markeredgecolor=pert_cors[g[0]][g[2]],
#                markeredgewidth=0.5, markerfacecolor='white', zorder=1000)
#     for _g, _d in pars[(pars['interval'] != 'median') & (pars['quantity'] == p)].groupby(['interval']):
#         ax[2].vlines(med_width['lower'], _d['lower'], _d['upper'])
plt.savefig('../../figures/Fig5_perturbation_protein.pdf',
            bbox_inches='tight')
plt.tight_layout()
plt.subplots_adjust(hspace=0.2)
