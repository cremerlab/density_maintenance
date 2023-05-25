# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Load datasets
growth_params = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
model_params = pd.read_csv(
    '../../data/mcmc/perturbation_parameter_percentiles.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ms_conf = pd.read_csv(
    '../../data/mass_spectrometry/periplasm_fractionation_confirmation.csv')
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
ppcs = ppcs[(ppcs['model'] == 'const_phi_mem') &
            (ppcs['volume_scale'] == 'linear_width')]
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']

# %%
err_widths = {'95%': 0.5, '75%': 1, '25%': 1.5}
fig, ax = plt.subplots(2, 2, figsize=(3, 2.25), sharex=True)
ax = ax.ravel()
for a in ax:
    a.set_xticks([0, 1, 2])
ax[0].set_ylim([0.3, 1.3])
ax[1].set_ylim([1, 6])
ax[2].set_ylim([-0.5, 5])
ax[3].set_ylim([0, 8])
ax[3].set_yticks([1, 3, 5, 7])
ax[1].set_yticks([1, 2,  3, 4, 5, 6])
ax[2].set_yticks([0, 1, 2, 3, 4])
ax[0].set_ylabel('average\nwidth [µm]', fontsize=6)
ax[1].set_ylabel('average\nlength [µm]', fontsize=6)
ax[2].set_ylabel('average\nvolume [µm$^3$]\n', fontsize=6)
ax[3].set_ylabel('average\naspect ratio [µm$^3$]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

pert_cors = {'wildtype': {'lacZ': cor['dark_black'],
                          'none': cor['primary_blue']},
             'malE-rbsB-fliC-KO': {'none': cor['primary_purple'],
                                   'malE': cor['primary_green'],
                                   'rbsB': cor['primary_gold']},
             'lpp14': {'none': cor['primary_red']}}

# Plot the percentiles
perc_cors = {'95%': cor['pale_black'], '75%': cor['light_black'],
             '25%': cor['primary_black'], 'median': cor['black']}
axes = {'w_rep': ax[0], 'ell_rep': ax[1],
        'vol_rep': ax[2], 'alpha_rep': ax[3]}
for i, (g, d) in enumerate(ppcs[ppcs['quantity'].isin(['w_rep', 'ell_rep',
                                                       'vol_rep', 'alpha_rep']) &
                                ppcs['interval'].isin(perc_cors.keys())
                                ].groupby(['quantity', 'interval'], sort=False)):
    if g[1] != 'median':
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], alpha=0.25,
                                color=perc_cors[g[1]], zorder=100)
    else:
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'], lw=1, alpha=0.25,
                        color=perc_cors[g[1]])

for g, d in size_data.groupby(['source']):
    for i, p in enumerate(['width_um', 'length_um', 'volume_um3', 'aspect_ratio']):
        ax[i].plot(d['growth_rate_hr'], d[p], mapper[g]['m'], ms=4,
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   color=mapper[g]['c'], alpha=0.45, zorder=500)

for g, d in growth_params.groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source', 'inducer', 'temperature'], sort=False):
    if (g[0] != 'wildtype') | (g[1] != 'none') | (g[2] != 0) | (g[5] != 37):
        continue
    sizes = model_params[(model_params['strain'] == g[0]) &
                         (model_params['overexpression'] == g[1]) &
                         (model_params['inducer_conc'] == g[2]) &
                         (model_params['carbon_source'] == g[3]) &
                         (model_params['quantity'].isin(['width_mu', 'length_mu', 'volume_mu', 'alpha_mu']))]
    med_growth = d[d['interval'] == 'median']
    for i, p in enumerate(['width_mu', 'length_mu', 'volume_mu', 'alpha_mu']):
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
plt.savefig('../../figures/Fig4_wildtype_dimensions.pdf', bbox_inches='tight')

# %%

fig, ax = plt.subplots(1, 2, figsize=(3.25, 1.75), sharex=True)
ax[0].set_ylim([0, 25])
ax[1].set_ylim([0, 0.125])
ax[0].set_xlim([0.2, 1.3])
ax[0].set_ylabel('$m_{peri}$ [fg / cell]', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$', fontsize=6)
ax[0].set_title('periplasmic protein mass', fontsize=6)
ax[1].set_title('allocation to periplasmic protein', fontsize=6)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)


for i, (g, d) in enumerate(ppcs[ppcs['quantity'].isin(['m_peri_rep', 'phi_peri_rep']) &
                                ppcs['interval'].isin(perc_cors.keys())
                                ].groupby(['quantity', 'interval'], sort=False)):
    if g[0] == 'm_peri_rep':
        a = ax[0]
    else:
        a = ax[1]
    if g[1] != 'median':
        a.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], alpha=0.25,
                       color=perc_cors[g[1]], zorder=100)
    else:
        a.plot(d['growth_rate_hr'], d['lower'], lw=1, alpha=0.25,
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

for g, d in growth_params.groupby(['strain', 'overexpression', 'inducer_conc', 'carbon_source'], sort=False):
    if (g[0] != 'wildtype') | (g[1] != 'none'):
        continue
    pars = model_params[(model_params['strain'] == g[0]) &
                        (model_params['overexpression'] == g[1]) &
                        (model_params['inducer_conc'] == g[2]) &
                        (model_params['carbon_source'] == g[3]) &
                        (model_params['quantity'].isin(['m_peri', 'phi_peri']))]
    if len(pars) == 0:
        continue
    med_growth = d[d['interval'] == 'median']
    for i, p in enumerate(['m_peri', 'phi_peri']):
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
for k, v in pert_cors.items():
    for _k, _v in v.items():
        ax[1].plot([], [], 'o', markeredgecolor=_v, markerfacecolor='w', markeredgewidth=1, ms=4,
                   label=f'{k}/{_k}')


plt.tight_layout()
plt.savefig('../../figures/Fig4_protein_protocol_validation.pdf',
            bbox_inches='tight')

# %%
# ######################################################################
# SUPPLEMENTAL FIGURES
# ######################################################################
fig, ax = plt.subplots(2, 1, figsize=(3, 3), sharex=True, sharey=True)

for a in ax:
    a.set_yticks([0, 1])
    a.set_yticklabels(['supernatant\nfraction', 'pellet\nfraction'])
    a.set_ylim([-0.5, 1.5])
axes = {'periplasmic': ax[0], 'non-periplasmic': ax[1]}

colors = ['purple', 'black']
ms_conf['sublocalization'] = 'periplasmic'
ms_conf.loc[ms_conf['localization'] != 'periplasm',
            'sublocalization'] = 'non-periplasmic'
for g, d in ms_conf.groupby(['sublocalization']):
    axes[g].set_title(
        f'annotated {g} proteins (N = {len(d)})', fontsize=6)
    print(g)
    _ = axes[g].plot(np.log2(d['relative_signal']), np.random.normal(
        0, 0.1, len(d)) + i, '.', color=cor[f'primary_{colors[i]}'], markeredgewidth=0, alpha=0.5,
        ms=3)

    # for i, p in enumerate(['peri_enrichment', 'cyto_enrichment']):
    # _ = axes[g].plot(np.log2(d[p]), np.random.normal(
    # 0, 0.1, len(d)) + i, '.', color=cor[f'primary_{colors[i]}'], markeredgewidth=0, alpha=0.5,
    # ms=3)

plt.tight_layout()
# plt.savefig('../../figures/FigSX_MS_periplasm_validation_plots.pdf')

# %%
# ######################################################################
# SUPPLEMENTAL FIGURES
# ######################################################################
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
