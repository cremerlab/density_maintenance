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

rbsB_oe = posts[(posts['strain'] == 'malE-rbsB-fliC-KO') &
                (posts['carbon_source'] == 'acetate') &
                (posts['overexpression'] == 'rbsB') &
                (posts['inducer_conc'] == 100)]

rbsB_oe_p = params[(params['strain'] == 'malE-rbsB-fliC-KO') &
                   (params['carbon_source'] == 'acetate') &
                   (params['overexpression'] == 'rbsB') &
                   (params['inducer_conc'] == 100)]

width_norm = malE_oe[malE_oe['parameter'] == 'width']['kde'].max()
aspect_norm = malE_oe[malE_oe['parameter'] == 'aspect_ratio']['kde'].max()
rho_norm = malE_oe[malE_oe['parameter'] == 'rho_ratio']['kde'].max()
norms = [width_norm, aspect_norm, rho_norm]
fig, ax = plt.subplots(1, 3, figsize=(3, 2), sharey=True)
ax[0].set_xlim([0.4, 0.9])
ax[1].set_xlim([2, 4])
err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}
# Plot the wildtypes
for i in reversed(range(4)):
    width = wt_acetate[wt_acetate['parameter'] == 'width_mu']
    aspect_ratio = wt_acetate[wt_acetate['parameter'] == 'aspect_ratio_mu']
    rho_ratio = wt_acetate[wt_acetate['parameter'] == 'rho_ratio']
    if i == 3:
        ls = '-'
    else:
        ls = '--'
    for j, prop in enumerate([width, aspect_ratio, rho_ratio]):
        if i == 3:
            ax[j].fill_between(prop['value'], 3 * np.ones(len(prop)), 3 + prop['kde'] / prop['kde'].max(), ls=ls,
                               lw=0.5, color=cor['pale_blue'], zorder=4 - i)
        ax[j].plot(prop['value'], i + prop['kde'] / prop['kde'].max(), ls=ls,
                   lw=0.5, color=cor['primary_blue'], zorder=4 - i)

for j in reversed(range(4)):
    for i, prop in enumerate(['width_mu', 'aspect_ratio_mu', 'rho_ratio']):
        _param = wt_acetate_p[wt_acetate_p['quantity'] == prop]
        _med = _param[_param['interval'] == 'median']
        _err = _param[_param['interval'] != 'median']
        for _g, _d in _err.groupby(['interval'], sort=False):
            ax[i].hlines(j + 0.5, _d['lower'], _d['upper'],
                         lw=err_widths[_g], color=cor['primary_blue'])
        ax[i].plot(_med['lower'], j + 0.5, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
                   markerfacecolor='w', zorder=1000)


# Plot the d3
for i, prop in enumerate(['width_mu', 'aspect_ratio_mu', 'rho_ratio']):
    q = d3_acetate[d3_acetate['parameter'] == prop]
    ax[i].plot(q['value'], 2 + q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['primary_purple'])
    ax[i].fill_between(q['value'], 2 * np.ones(len(q)),  2 +
                       q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_purple'], alpha=0.5)

    _param = d3_acetate_p[d3_acetate_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(2.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['primary_purple'])
    ax[i].plot(_med['lower'], 2.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_purple'],
               markerfacecolor='w', zorder=1000)

for i, prop in enumerate(['width_mu', 'aspect_ratio_mu', 'rho_ratio']):
    q = malE_oe[malE_oe['parameter'] == prop]
    ax[i].plot(q['value'], q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['green'])
    ax[i].fill_between(q['value'], np.zeros(len(q)),
                       q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_green'], alpha=0.5)
    _param = malE_oe_p[malE_oe_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(0.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['green'])
    ax[i].plot(_med['lower'], 0.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['green'],
               markerfacecolor='w', zorder=1000)


for i, prop in enumerate(['width_mu', 'aspect_ratio_mu', 'rho_ratio']):
    q = rbsB_oe[rbsB_oe['parameter'] == prop]
    ax[i].plot(q['value'], 1 + q['kde'] / q['kde'].max(), '-',
               lw=0.5, color=cor['gold'])
    ax[i].fill_between(q['value'],  np.ones(len(q)),
                       1 + q['kde'] / q['kde'].max(), lw=0.5, color=cor['pale_gold'], alpha=0.5)
    _param = rbsB_oe_p[rbsB_oe_p['quantity'] == prop]
    _med = _param[_param['interval'] == 'median']
    _err = _param[_param['interval'] != 'median']
    for _g, _d in _err.groupby(['interval'], sort=False):
        ax[i].hlines(1.4, _d['lower'], _d['upper'],
                     lw=err_widths[_g], color=cor['gold'])
    ax[i].plot(_med['lower'], 1.4, 'o', ms=3, markeredgewidth=0.5, markeredgecolor=cor['gold'],
               markerfacecolor='w', zorder=1000)

for a in ax:
    a.set_yticks([])
plt.savefig('../../figures/Fig3_periplasm_oe.pdf', bbox_inches='tight')
