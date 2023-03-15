# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import Blues_7
import size.viz
import size.analytical
cor, pal = size.viz.matplotlib_style()

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')

# Restrict to wildtype
params = params[~((params['overexpression'] != 'none')
                  * (params['overexpression'] == 0))]
posts = posts[~((posts['overexpression'] != 'none')
                * (posts['overexpression'] == 0))]
params.loc[params['quantity'] == 'prot_per_cell', 'lower'] *= (1E-6/1E-15)
params.loc[params['quantity'] == 'prot_per_cell', 'upper'] *= (1E-6/1E-15)
posts.loc[posts['parameter'] == 'prot_per_cell', 'value'] *= (1E-6/1E-15)
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'}, inplace=True)
shapes = shapes[~((shapes['overexpression'] != 'none')
                  * (shapes['overexpression'] == 0))]
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']

# %%

# Separate the samples by hand so I know what I'm doing
wt_acetate = posts[(posts['strain'] == 'wildtype') &
                   (posts['carbon_source'] == 'acetate') &
                   (posts['overexpression'] == 'none') &
                   (posts['inducer_conc'] == 0)]

wt_acetate_p = params[(params['strain'] == 'wildtype') &
                      (params['carbon_source'] == 'acetate') &
                      (params['overexpression'] == 'none') &
                      (params['inducer_conc'] == 0)]

d3_acetate = posts[(posts['strain'] == 'malE-rbsB-fliC-KO') &
                   (posts['carbon_source'] == 'acetate') &
                   (posts['overexpression'] == 'none') &
                   (posts['inducer_conc'] == 0)]
d3_acetate_p = params[(params['strain'] == 'malE-rbsB-fliC-KO') &
                      (params['carbon_source'] == 'acetate') &
                      (params['overexpression'] == 'none') &
                      (params['inducer_conc'] == 0)]


malE_oe = posts[(posts['strain'] == 'malE-rbsB-fliC-KO') &
                (posts['carbon_source'] == 'acetate') &
                (posts['overexpression'] == 'malE') &
                (posts['inducer_conc'] == 100)]

malE_oe_p = params[(params['strain'] == 'malE-rbsB-fliC-KO') &
                   (params['carbon_source'] == 'acetate') &
                   (params['overexpression'] == 'malE') &
                   (params['inducer_conc'] == 100)]

# rbsB_oe = posts[(posts['strain'] == 'malE-rbsB-fliC-KO') &
#                 (posts['carbon_source'] == 'acetate') &
#                 (posts['overexpression'] == 'rbsB') &
#                 (posts['inducer_conc'] == 100)]

# rbsB_oe_p = params[(params['strain'] == 'malE-rbsB-fliC-KO') &
#                    (params['carbon_source'] == 'acetate') &
#                    (params['overexpression'] == 'rbsB') &
#                    (params['inducer_conc'] == 100)]

# width_norm = malE_oe[malE_oe['parameter'] == 'width']['kde'].max()
# aspect_norm = malE_oe[malE_oe['parameter'] == 'aspect_ratio']['kde'].max()
# rho_norm = malE_oe[malE_oe['parameter'] == 'rho_ratio']['kde'].max()
# norms = [width_norm, aspect_norm, rho_norm]
fig, ax = plt.subplots(1, 3, figsize=(3, 2), sharey=True)
ax[0].set_xlim([0.02, 0.06])
ax[1].set_xlim([0.45, 0.85])
ax[2].set_xlim([2, 4])
ax[0].set_xticks([0.02, 0.04, 0.06])
ax[1].set_xticks([0.5, 0.6, 0.7, 0.8])
ax[2].set_xticks([2, 2.5, 3, 3.5, 4])
for a in ax:
    a.xaxis.tick_top()

err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}
# Plot the wildtypes
for i in reversed(range(3)):
    width = wt_acetate[wt_acetate['parameter'] == 'width_mu']
    aspect_ratio = wt_acetate[wt_acetate['parameter'] == 'aspect_ratio_mu']
    phi_M = wt_acetate[wt_acetate['parameter'] == 'phi_M']
    if i == 2:
        ls = '-'
        c = cor['primary_blue']
    else:
        ls = '--'
        c = cor['light_blue']
    for j, prop in enumerate([phi_M, width, aspect_ratio]):
        if i == 2:
            ax[j].fill_between(prop['value'], (i+1) * np.ones(len(prop)), (i + 1) + prop['kde'] / prop['kde'].max(), ls=ls,
                               lw=0.5, color=cor['pale_blue'], zorder=4 - i)
        # if i >= 1:
        ax[j].plot(prop['value'], (i+1) + prop['kde'] / prop['kde'].max(), ls=ls,
                   lw=0.5, color=c, zorder=4 - i)


for i, prop in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
    _param = wt_acetate_p[wt_acetate_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        for j in range(3):
            ax[i].hlines(j + 1.5, _d['lower'], _d['upper'],
                         lw=err_widths[_g], color=cor['primary_blue'])
            ax[i].plot(_med['lower'], j + 1.5, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
                       markerfacecolor='w', zorder=1000)


# Plot the d3
for i, prop in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
    q = d3_acetate[d3_acetate['parameter'] == prop]
    ax[i].plot(q['value'], 2 + q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['primary_purple'])
    ax[i].fill_between(q['value'], 2 * np.ones(len(q)),  2 +
                       q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_purple'], alpha=0.5)
    # ax[i].plot(q['value'], 1 + q['kde'] / q['kde'].max(), '--',
    #    lw=0.5, color=cor['primary_purple'])

    _param = d3_acetate_p[d3_acetate_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(2.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['primary_purple'])
    ax[i].plot(_med['lower'], 2.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_purple'],
               markerfacecolor='w', zorder=1000)

for i, prop in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
    q = malE_oe[malE_oe['parameter'] == prop]
    ax[i].plot(q['value'], 1 + q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['green'])
    ax[i].fill_between(q['value'], np.ones(len(q)),
                       1 + q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_green'], alpha=0.5)
    _param = malE_oe_p[malE_oe_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(1.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['green'])
    ax[i].plot(_med['lower'], 1.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['green'],
               markerfacecolor='w', zorder=1000)


# for i, prop in enumerate(['width_mu', 'aspect_ratio_mu', 'rho_ratio']):
#     q = rbsB_oe[rbsB_oe['parameter'] == prop]
#     ax[i].plot(q['value'], 1 + q['kde'] / q['kde'].max(), '-',
#                lw=0.5, color=cor['gold'])
#     ax[i].fill_between(q['value'],  np.ones(len(q)),
#                        1 + q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_gold'], alpha=0.5)
#     _param = rbsB_oe_p[rbsB_oe_p['quantity'] == prop]
#     _med = _param[_param['interval'] == 'median']
#     _err = _param[_param['interval'] != 'median']
#     for _g, _d in _err.groupby(['interval'], sort=False):
#         ax[i].hlines(1.4, _d['lower'], _d['upper'],
#                      lw=err_widths[_g], color=cor['gold'])
#     ax[i].plot(_med['lower'], 1.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['gold'],
#                markerfacecolor='w', zorder=1000)

for a in ax:
    a.set_yticks([])
plt.savefig('../../figures/Fig3_periplasm_oe.pdf', bbox_inches='tight')

# %%
wt_glucose = posts[(posts['strain'] == 'wildtype') &
                   (posts['carbon_source'] == 'glucose') &
                   (posts['overexpression'] == 'none') &
                   (posts['inducer_conc'] == 0)]

wt_glucose_p = params[(params['strain'] == 'wildtype') &
                      (params['carbon_source'] == 'glucose') &
                      (params['overexpression'] == 'none') &
                      (params['inducer_conc'] == 0)]

lacZ_glucose = posts[(posts['strain'] == 'wildtype') &
                     (posts['carbon_source'] == 'glucose') &
                     (posts['overexpression'] == 'lacZ') &
                     (posts['inducer_conc'] == 20)]

lacZ_glucose_p = params[(params['strain'] == 'wildtype') &
                        (params['carbon_source'] == 'glucose') &
                        (params['overexpression'] == 'lacZ') &
                        (params['inducer_conc'] == 20)]

fig, ax = plt.subplots(1, 3, figsize=(3, 1.3), sharey=True)
ax[0].set_xlim([0.005, 0.035])
ax[1].set_xlim([0.65, 0.95])
ax[2].set_xlim([2.5, 4.5])
# ax[0].set_xticks([0, 0.02, 0.04])
ax[1].set_xticks([0.7, 0.8, 0.9])
ax[2].set_xticks([3, 3.5, 4])
for i in reversed(range(2)):
    width = wt_glucose[wt_glucose['parameter'] == 'width_mu']
    aspect_ratio = wt_glucose[wt_glucose['parameter'] == 'aspect_ratio_mu']
    phi_M = wt_glucose[wt_glucose['parameter'] == 'phi_M']
    if i == 2:
        ls = '-'
        c = cor['primary_blue']
    else:
        ls = '--'
        c = cor['light_blue']
    for j, prop in enumerate([phi_M, width, aspect_ratio]):
        if i == 1:
            ax[j].fill_between(prop['value'], i * np.ones(len(prop)), i + prop['kde'] / prop['kde'].max(), ls=ls,
                               lw=0.5, color=cor['pale_blue'], zorder=4 - i)
        ax[j].plot(prop['value'], i + prop['kde'] / prop['kde'].max(), ls=ls,
                   lw=0.5, color=c, zorder=4 - i)


for i, prop in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
    _param = wt_glucose_p[wt_glucose_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        for j in range(2):
            ax[i].hlines(j + 0.5, _d['lower'], _d['upper'],
                         lw=err_widths[_g], color=cor['primary_blue'], zorder=1000)
            ax[i].plot(_med['lower'], j + 0.5, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
                       markerfacecolor='w', zorder=1005)

for i, prop in enumerate(['phi_M', 'width_mu', 'aspect_ratio_mu']):
    q = lacZ_glucose[lacZ_glucose['parameter'] == prop]
    ax[i].plot(q['value'], q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['primary_black'])
    ax[i].fill_between(q['value'], np.zeros(len(q)),
                       q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_black'], alpha=0.5)
    _param = lacZ_glucose_p[lacZ_glucose_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(0.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['primary_black'])
    ax[i].plot(_med['lower'], 0.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
               markerfacecolor='w', zorder=1000)

for a in ax:
    a.set_yticks([])
    # a.xaxis.tick_top()
plt.savefig('../../figures/Fig3_cytoplasm_oe.pdf', bbox_inches='tight')

# %%
lacZ_glucose_p
