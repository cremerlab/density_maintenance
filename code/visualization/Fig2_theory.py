# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.analytical
import seaborn as sns
mapper = size.viz.lit_mapper()
err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}
cor, pal = size.viz.matplotlib_style()
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
# mass_fracs = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')
pred_phi = pd.read_csv('../../data/mcmc/predicted_phi_scaling_lppwt.csv')
singular_posts = pd.read_csv('../../data/mcmc/singular_posterior_kde.csv')
singular_percs = pd.read_csv(
    '../../data/mcmc/singular_parameter_percentiles.csv')
mass_spec_medians = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')
mass_spec_medians = mass_spec_medians[mass_spec_medians['interval'] == 'median']

# Restrict to wildtype
params = params[(params['strain'] == 'wildtype') &
                (params['overexpression'] == 'none') &
                (params['inducer_conc'] == 0)]
posts = posts[(posts['strain'] == 'wildtype') &
              (posts['overexpression'] == 'none') &
              (posts['inducer_conc_ng_mL'] == 0)]
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'})
shapes = shapes[(shapes['strain'] == 'wildtype') &
                (shapes['overexpression'] == 'none') &
                (shapes['inducer_conc'] == 0)]
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']


# %%
# Define the parameters
k = np.array([1, 2, 3])
delta = 0.025
slope = k * delta
phi_range = np.linspace(0.01, 0.10, 100)
sav_range = np.linspace(5.8001, 7.2, 100)

# %%

fig, ax = plt.subplots(1, 1, figsize=(2, 1.75))

ls = ['-', '--', ':', '-.']
for i, s in enumerate(slope):
    theo = s * (sav_range - 5.8)
    ax.plot(sav_range, theo + 0.02,  ls=ls[i],
            lw=1, color=cor['primary_blue'], zorder=1000)

# Plot the mass spec data
for g, d in mass_spec_medians.groupby(['dataset_name']):
    ax.plot(d[d['quantity'] == 'mass_spec_sav']['lower'],
            d[d['quantity'] == 'mass_spec_phi_M']['lower'],
            linestyle='none', marker=mapper[g]['m'], markeredgecolor='k',
            alpha=0.15, markeredgewidth=0.5, color=mapper[g]['c'], ms=4, zorder=10)


med_params = params[(params['quantity'].isin(['surface_area_vol_mu', 'phi_M'])) &
                    (params['interval'] == 'median')]
for g, d in params[params['quantity'].isin(['surface_area_vol_mu', 'phi_M'])].groupby(['carbon_source',
                                                                                       'interval'], sort=False):
    if g[0] == 'LB':
        continue

    if g[1] == 'median':
        ax.plot(d[d['quantity'] == 'surface_area_vol_mu']['lower'], d[d['quantity'] == 'phi_M']['upper'],
                marker='o', markeredgewidth=0.75, markeredgecolor=cor['primary_blue'],
                markerfacecolor='white', ms=3, zorder=1000)
    else:
        sav = d[d['quantity'] == 'surface_area_vol_mu']
        phiM = d[d['quantity'] == 'phi_M']
        med_sav = med_params[(med_params['quantity'] == 'surface_area_vol_mu') &
                             (med_params['carbon_source'] == g[0])]
        med_phi_M = med_params[(med_params['quantity'] == 'phi_M') &
                               (med_params['carbon_source'] == g[0])]
        ax.vlines(med_sav['lower'], phiM['lower'], phiM['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)
        ax.hlines(med_phi_M['lower'], sav['lower'], sav['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)

ax.set_xlim([5.5, 7.2])
# ax.set_xticks([5, 5.5, 6, 6.5, 7])
ax.set_ylim([0.01, 0.12])
ax.set_xlabel('surface area to volume [µm$^{-1}$]', fontsize=6)
ax.set_ylabel(
    'periplasmic biomass fraction', fontsize=6)
plt.savefig('../../figures/Fig2_sav_scaling.pdf', bbox_inches='tight')


# %%
# %%
# Set up the canvas and label
fig, ax = plt.subplots(1, 3, figsize=(2, 1.5))
for a in ax:
    a.set_yticks([])
ax[0].set_xlim([0.5, 0.9])
ax[1].set_xlim([1.8, 3.2])
ax[2].set_xlim([2.5, 4])
nudge = 1

axes = {'width_mu': 0, 'length_mu': 1, 'aspect_ratio_mu': 2}
locs = {'glucoseCAA': 5, 'glucose': 4,
        'glycerol': 3, 'sorbitol': 2, 'acetate': 1}
for i, (g, d) in enumerate(posts[posts['parameter'].isin(['width_mu', 'length_mu', 'aspect_ratio_mu'])].groupby(['carbon_source'])):
    if g == 'LB':
        continue
    for _g, _d in d.groupby(['parameter']):
        ax[axes[_g]].fill_between(
            _d['value'], locs[g] * nudge * np.ones(len(_d)),
            locs[g] * nudge * np.ones(len(_d)) + _d['kde'] / _d['kde'].max(),
            color=cor['pale_blue'], alpha=0.4, zorder=i+1)
        ax[axes[_g]].plot(_d['value'], locs[g] * nudge + _d['kde'] / _d['kde'].max(), '-',
                          color=cor['primary_blue'], lw=0.75, zorder=i+1)

for g, d in singular_posts[singular_posts['parameter'].isin(['alpha'])].groupby(['parameter']):
    ax[axes[g]].fill_between(d['value'], 0, d['kde'] /
                             d['kde'].max(), color=cor['pale_black'])
    ax[axes[g]].plot(d['value'], d['kde'] / d['kde'].max(),
                     color=cor['primary_black'], lw=0.75)

err_widths = {'95%': 0.25, '75%': 1, '25%': 2.5}

for g, d in errs[errs['quantity'].isin(list(axes.keys()))
                 ].groupby(['carbon_source', 'quantity', 'interval']):
    if g[0] == 'LB':
        continue
    ax[axes[g[1]]].hlines(locs[g[0]] * nudge + 0.5, d['lower'], d['upper'], lw=err_widths[g[-1]],
                          color=cor['primary_blue'], zorder=1000)


for g, d in medians[medians['quantity'].isin(list(axes.keys()))
                    ].groupby(['carbon_source', 'quantity']):
    if g[0] == 'LB':
        continue
    ax[axes[g[1]]].plot(d['lower'], locs[g[0]] * nudge + 0.5, 'o', markeredgewidth=0.5, ms=2,
                        markeredgecolor=cor['primary_blue'], markerfacecolor='white',
                        zorder=1000)

for g, d in singular_errs[singular_errs['quantity'].isin(list(axes.keys()))
                          ].groupby(['quantity', 'interval']):
    ax[axes[g[0]]].hlines(0.5, d['lower'], d['upper'], lw=err_widths[g[-1]],
                          color=cor['primary_black'], zorder=1000)

for g, d in singular_medians[singular_medians['quantity'].isin(list(axes.keys()))
                             ].groupby(['quantity', 'interval']):
    ax[axes[g[0]]].plot(d['lower'], 0.5, 'o', markeredgewidth=0.5,
                        markerfacecolor='w', markersize=3, markeredgecolor=cor['primary_black'], zorder=1000)

plt.tight_layout()
plt.savefig('../../figures/Fig2_constant_parameters_kde.pdf')
# %%
# mass_fracs = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
# mass_fracs = mass_fracs[mass_fracs['periplasm']]
# dry_frac = 0.3
# prot_frac = 0.65
# density = 1.1

# # %%
# # Do the proper classification
# # genes = pd.read_csv('../../data/literature/genes_classification_all.csv')
# # _genes = genes[genes['location'].isin(['IM', 'OM', 'PE', 'LPO'])]
# # mass_fracs = mass_fracs[mass_fracs['gene_name'].isin(_genes['gene'].unique())]
lit_size = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
# growth = pd.read_csv('../../data/summaries/summarized_growth_measurements.csv')
# growth = growth[(growth['strain'] == 'wildtype') & (
#     growth['overexpression'] == 'none') & (growth['inducer_conc'] == 0)]
# growth = growth.groupby(['carbon_source']).mean().reset_index()
# wt_size = size_data[(size_data['strain'] == 'wildtype') & (
#     size_data['overexpression'] == 'none') & (size_data['inducer_conc'] == 0)]
# for g, d in growth.groupby(['carbon_source', 'growth_rate_hr']):
#     wt_size.loc[wt_size['carbon_source'] == g[0], 'growth_mu'] = g[1]

# wt_size
# # %%
# lit_size.dropna(inplace=True)
# wt_size.dropna(inplace=True)
# # growth = np.concatenate(
# # [lit_size['growth_rate_hr'].values, wt_size['growth_mu'].values]).flatten()
# growth = wt_size['growth_mu']
# # ell = np.concatenate([lit_size['length_um'].values,
# #  wt_size['length'].values]).flatten()
# ell = wt_size['length']
# # w = np.concatenate([lit_size['width_um'].values,
# #    wt_size['width_median'].values]).flatten()
# w = wt_size['width_median']
# # v = np.concatenate(
# # [lit_size['volume_um3'], wt_size['volume'].values]).flatten()
# v = wt_size['volume']
# # sav = np.concatenate([lit_size['surface_to_volume'].values,
# #  wt_size['surface_to_volume'].values]).flatten()
# sav = wt_size['surface_to_volume']
# peri_vol = np.pi * ell * w * 0.025


# %%
# # Determine the simple relations
# w_popt = scipy.stats.linregress(growth, w)
# ell_popt = scipy.stats.linregress(growth, ell)
# peri_vol_popt = scipy.stats.linregress(growth, peri_vol)
# vol_popt = scipy.stats.linregress(growth, v)
# sav_popt = scipy.stats.linregress(growth, sav)

# # Compute the periplasmic protein density
# mass_fracs['width'] = w_popt[0] * \
#     mass_fracs['growth_rate_hr'] + w_popt[1]
# mass_fracs['length'] = np.exp(
#     ell_popt[0] * mass_fracs['growth_rate_hr'] + ell_popt[1])
# mass_fracs['peri_vol'] = np.exp(
#     peri_vol_popt[0] * mass_fracs['growth_rate_hr'] + peri_vol_popt[1])
# mass_fracs['volume'] = np.exp(
#     vol_popt[0] * mass_fracs['growth_rate_hr'] + vol_popt[1])
# mass_fracs['sav'] = np.exp(
#     sav_popt[0] * mass_fracs['growth_rate_hr'] + sav_popt[1])
# # mass_fracs['peri_volume'] = size.analytical.surface_area(mass_fracs['length'], mass_fracs['width']) * 0.025
# mass_fracs['tot_protein'] = density * \
#     dry_frac * prot_frac * mass_fracs['volume']
# mass_fracs['peri_protein'] = mass_fracs['mass_frac'] * \
#     mass_fracs['tot_protein']
# mass_fracs['rho_peri'] = (
#     mass_fracs['peri_protein'] * 1E3) / mass_fracs['peri_vol']
# mass_fracs['biomass_frac'] = mass_fracs['peri_protein'] / \
#     (density * dry_frac * mass_fracs['volume'])
# mass_fracs = mass_fracs.groupby(
#     ['dataset_name', 'condition', 'growth_rate_hr', 'width']).sum().reset_index()

markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h',
           'p', 'P', '*', 'o', '8', 'd', '>', 'v', '<', '^']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(lit_size['source'].unique())
for n in mass_spec_medians['dataset_name'].unique():
    names.append(n)
mapper = {n: {'m': m, 'c': c} for n, m, c in zip(names, markers, cors)}
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.75), sharey=True)

# for g, d in mass_fracs.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
#     ax.plot(d['biomass_frac'], d['width'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
#             markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
#             alpha=0.75, ms=3)


# Plot the percentiles and the medians
for g, d in errs[(errs['quantity'].isin(['phi_M', 'width_mu']))
                 ].groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    phi = medians[(medians['carbon_source'] == g[0]) &
                  (medians['quantity'] == 'phi_M')]['lower']
    w = medians[(medians['carbon_source'] == g[0]) & (
        medians['quantity'] == 'width_mu')]['lower']
    phi_d = d[d['quantity'] == 'phi_M']
    w_d = d[d['quantity'] == 'width_mu']
    ax.hlines(w, phi_d['lower'], phi_d['upper'], lw=err_widths[g[1]],
              color=cor['primary_blue'])
    ax.vlines(phi, w_d['lower'], w_d['upper'], lw=err_widths[g[1]],
              color=cor['primary_blue'])

for g, d in medians[medians['quantity'].isin(['phi_M',
                                              'width_mu'])].groupby(['carbon_source']):
    if g == 'LB':
        continue
    phi = d[d['quantity'] == 'phi_M']['lower']
    w = d[d['quantity'] == 'width_mu']['lower']
    ax.plot(phi, w, 'o', markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
            markerfacecolor='white', ms=3)


band_cor = {'95%': '#EACFCE',
            '75%': '#E0B1B1',
            '25%': '#D38888'}
for i, (g, d) in enumerate(pred_err[pred_err['quantity'] == 'width_simple'].groupby(['interval'], sort=False)):
    ax.fill_between(d['phi_M'], d['lower'], d['upper'],
                    color=band_cor[g])
ax.set_ylim([0.45, 1])
ax.set_xlim([0, 0.06])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_ylabel('w\naverage width [µm]', fontsize=6)
plt.savefig('../../figures/Fig2_width_prediction_wildtype.pdf')
# %%
phi_range = np.linspace(0, 0.06, 10)
k = 0.17
delta = 0.024
theo = (1 + (k / (phi_range**-1 - 1)))**-1
plt.plot(phi_range, theo)
# plt.ylim([0, 10])

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))

sav_range = np.linspace(4, 10, 100)
x = [0.1, 0.2, 0.45]
delta = 0.025
for i, _x in enumerate(x):
    if i == 0:
        ls = '--'
    elif i == 1:
        ls = '-.'
    else:
        ls = ':'
    theo = (_x**-1 * ((sav_range * delta)**-1 - 1) + 1)**-1
    ax.plot(sav_range, theo, ls=ls, lw=1, color=cor['light_black'])


for g, d in mass_spec_medians.groupby(['dataset_name']):
    ax.plot(d[d['quantity'] == 'mass_spec_sav']['lower'], d[d['quantity'] == 'mass_spec_phi_M']['lower'],
            label=g, marker=mapper[g]['m'], color=mapper[g]['c'], linestyle='none',
            ms=3, markeredgecolor=cor['primary_black'], markeredgewidth=0.25)

for g, d in errs.groupby(['carbon_source', 'interval']):
    if g[0] == 'LB':
        continue
    phi_loc = medians[(medians['quantity'] == 'phi_M') & (
        medians['carbon_source'] == g[0])]['lower']
    sv_loc = medians[(medians['quantity'] == 'surface_area_vol_mu') &
                     (medians['carbon_source'] == g[0])]['lower']
    phi = d[d['quantity'] == 'phi_M']
    sv = d[d['quantity'] == 'surface_area_vol_mu']
    ax.vlines(sv_loc, phi['lower'], phi['upper'], color=cor['primary_blue'],
              lw=err_widths[g[1]])
    ax.hlines(phi_loc, sv['lower'], sv['upper'], color=cor['primary_blue'],
              lw=err_widths[g[1]])

for g, d in medians.groupby(['carbon_source']):
    if g == 'LB':
        continue
    ax.plot(d[d['quantity'] == 'surface_area_vol_mu']['lower'],
            d[d['quantity'] == 'phi_M']['lower'], 'o',
            ms=3, markeredgewidth=0.5, markeredgecolor=cor['primary_blue'],
            markerfacecolor='w')

ax.set_ylim([0, 0.06])
ax.set_ylabel('periplasmic biomass fraction\n$\phi_{M}$', fontsize=6)
ax.set_xlabel('$S/V$\nsurface area to volume [µm$^{-1}$]', fontsize=6)
plt.savefig('../../figures/Fig2_SV_scaling.pdf', bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
width_range = np.linspace(0.45, 1, 100)
pref = (12 * 3 * 3 * 0.025) / (3 * 3 - 1)
w_min = 0.6
phi_max = 0.08
beta_0 = phi_max + pref / w_min
theo = phi_max - pref * (width_range - w_min)


for i, _k in enumerate([1, 2, 3]):
    pref = (12 * 3 * _k * 0.025) / (3 * 3 - 1)
    w_min = 0.61
    phi_max = 0.08
    theo = phi_max - pref * (width_range - w_min)

    ax.plot(width_range, theo, ls=ls[i],
            lw=1, color=cor['primary_blue'], zorder=100)

# Plot the mass spec data
for g, d in mass_spec_medians.groupby(['dataset_name']):
    ax.plot(d[d['quantity'] == 'mass_spec_widths']['lower'],
            d[d['quantity'] == 'mass_spec_phi_M']['lower'],
            linestyle='none', marker=mapper[g]['m'], markeredgecolor='k',
            alpha=0.35, markeredgewidth=0.5, color=mapper[g]['c'], ms=4, zorder=10)


med_params = params[(params['quantity'].isin(['width_mu', 'phi_M'])) &
                    (params['interval'] == 'median')]
for g, d in params[params['quantity'].isin(['width_mu', 'phi_M'])].groupby(['carbon_source',
                                                                            'interval'], sort=False):
    if g[0] == 'LB':
        continue

    if g[1] == 'median':
        ax.plot(d[d['quantity'] == 'width_mu']['lower'], d[d['quantity'] == 'phi_M']['upper'],
                marker='o', markeredgewidth=0.75, markeredgecolor=cor['primary_blue'],
                markerfacecolor='white', ms=3, zorder=1000)
    else:
        sav = d[d['quantity'] == 'width_mu']
        phiM = d[d['quantity'] == 'phi_M']
        med_sav = med_params[(med_params['quantity'] == 'surface_area_vol_mu') &
                             (med_params['carbon_source'] == g[0])]
        med_phi_M = med_params[(med_params['quantity'] == 'phi_M') &
                               (med_params['carbon_source'] == g[0])]
        ax.vlines(med_sav['lower'], phiM['lower'], phiM['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)
        ax.hlines(med_phi_M['lower'], sav['lower'], sav['upper'], lw=err_widths[g[1]],
                  color=cor['primary_blue'], zorder=99)

ax.set_xlim([0.6, 1])
# ax.set_xticks([5, 5.5, 6, 6.5, 7])
ax.set_ylim([0, 0.1])
ax.set_xlabel('average cell width [µm]', fontsize=6)
ax.set_ylabel('periplasmic protein biomass fraction', fontsize=6)
plt.savefig('../../figures/Fig2_width_scaling.pdf', bbox_inches='tight')


# %%

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')
singular_posts = pd.read_csv('../../data/mcmc/singular_posterior_kde.csv')
params = params[~((params['overexpression'] != 'none') &
                (params['inducer_conc'] == 0))]
posts = posts[~((posts['overexpression'] != 'none') &
              (posts['inducer_conc_ng_mL'] == 0))]
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'})
shapes = shapes[~((shapes['overexpression'] != 'none') &
                (shapes['inducer_conc'] == 0))]
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']


fig, ax = plt.subplots(1, 1, figsize=(2, 1.75), sharey=True)

cond_colors = {('wildtype', 'none'): cor['primary_blue'],
               ('malE-rbsB-fliC-KO', 'none'): cor['primary_purple'],
               ('malE-rbsB-fliC-KO', 'malE'): cor['primary_green'],
               ('wildtype', 'lacZ'): cor['primary_black']}

# Plot the percentiles and the medians
for g, d in errs[(errs['quantity'].isin(['phi_M', 'width_mu']))
                 ].groupby(['carbon_source', 'interval', 'strain', 'overexpression', 'inducer_conc']):

    if (g[0] == 'LB'):
        continue
    c = cond_colors[(g[2], g[3])]
    phi = medians[(medians['carbon_source'] == g[0]) &
                  (medians['strain'] == g[2]) &
                  (medians['overexpression'] == g[3]) &
                  (medians['inducer_conc'] == g[4]) &
                  (medians['quantity'] == 'phi_M')]['lower']
    w = medians[(medians['carbon_source'] == g[0]) &
                (medians['quantity'] == 'width_mu') &
                (medians['strain'] == g[2]) &
                (medians['overexpression'] == g[3]) &
                (medians['inducer_conc'] == g[4])]['lower']
    if (len(phi) >= 1) & (len(w) >= 1):
        phi_d = d[d['quantity'] == 'phi_M']
        w_d = d[d['quantity'] == 'width_mu']
        ax.hlines(w, phi_d['lower'], phi_d['upper'], lw=err_widths[g[1]],
                  color=c)
        ax.vlines(phi, w_d['lower'], w_d['upper'], lw=err_widths[g[1]],
                  color=c)

for g, d in medians[medians['quantity'].isin(['phi_M',
                                              'width_mu'])].groupby(['carbon_source', 'strain', 'overexpression', 'inducer_conc']):
    if (g[0] == 'LB') | (g[2] == 'rbsB'):
        continue
    c = cond_colors[(g[1], g[2])]
    phi = d[d['quantity'] == 'phi_M']['lower']
    w = d[d['quantity'] == 'width_mu']['lower']
    if (len(phi) >= 1) & (len(w) >= 1):
        ax.plot(phi, w, 'o', markeredgewidth=0.5, markeredgecolor=c,
                markerfacecolor='white', ms=3)

for g, d in mass_fracs.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    ax.plot(d['biomass_frac'], d['width'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
            markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
            alpha=0.5, ms=3)

band_cor = {'95%': '#EACFCE',
            '75%': '#E0B1B1',
            '25%': '#D38888'}
w_range = np.linspace(0.5, 1.2)
k = 0.1
alpha = 3
delta = 0.025
Lam = 12 * alpha * delta / (3 * alpha - 1)
theo = ((w_range - k) / (Lam * k) + 1)**-1
# for i, (g, d) in enumerate(pred_err[pred_err['quantity'] == 'width_simple'].groupby(['interval'], sort=False)):
# ax.fill_between(d['phi_M'], d['lower'], d['upper'],
# color=band_cor[g])
plt.plot(w_range, theo, '-')
ax.set_ylim([0.45, 1])
ax.set_xlim([0, 0.06])
ax.set_xlabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_ylabel('w\naverage width [µm]', fontsize=6)
plt.savefig('../../figures/Fig2_width_prediction_oe.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 1.75), sharey=True)

cond_colors = {('wildtype', 'none'): cor['primary_blue'],
               ('malE-rbsB-fliC-KO', 'none'): cor['primary_purple'],
               ('malE-rbsB-fliC-KO', 'malE'): cor['primary_green'],
               ('malE-rbsB-fliC-KO', 'rbsB'): cor['primary_gold'],
               ('wildtype', 'lacZ'): cor['primary_black']}

# Plot the percentiles and the medians
for g, d in errs[(errs['quantity'].isin(['phi_M', 'width_mu']))
                 ].groupby(['carbon_source', 'interval', 'strain', 'overexpression', 'inducer_conc']):

    if (g[0] == 'LB'):  # | (g[3] == 'rbsB'):
        continue
    c = cond_colors[(g[2], g[3])]
    phi = medians[(medians['carbon_source'] == g[0]) &
                  (medians['strain'] == g[2]) &
                  (medians['overexpression'] == g[3]) &
                  (medians['inducer_conc'] == g[4]) &
                  (medians['quantity'] == 'phi_M')]['lower']
    w = medians[(medians['carbon_source'] == g[0]) &
                (medians['quantity'] == 'width_mu') &
                (medians['strain'] == g[2]) &
                (medians['overexpression'] == g[3]) &
                (medians['inducer_conc'] == g[4])]['lower']
    if (len(phi) >= 1) & (len(w) >= 1):
        phi_d = d[d['quantity'] == 'phi_M']
        w_d = d[d['quantity'] == 'width_mu']
        ax.vlines(w, phi_d['lower'], phi_d['upper'], lw=err_widths[g[1]],
                  color=c)
        ax.hlines(phi, w_d['lower'], w_d['upper'], lw=err_widths[g[1]],
                  color=c)

for g, d in medians[medians['quantity'].isin(['phi_M',
                                              'width_mu'])].groupby(['carbon_source', 'strain', 'overexpression', 'inducer_conc']):
    if (g[0] == 'LB'):  # | (g[2] == 'rbsB'):
        continue
    c = cond_colors[(g[1], g[2])]
    phi = d[d['quantity'] == 'phi_M']['lower']
    w = d[d['quantity'] == 'width_mu']['lower']
    if (len(phi) >= 1) & (len(w) >= 1):
        ax.plot(w, phi, 'o', markeredgewidth=0.5, markeredgecolor=c,
                markerfacecolor='white', ms=3)

for g, d in mass_fracs.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    ax.plot(d['width'], d['biomass_frac'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
            markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
            alpha=0.5, ms=3, zorder=1)

band_cor = {'95%': '#EACFCE',
            '75%': '#E0B1B1',
            '25%': '#D38888'}
# w_range = np.linspace(0.5, 1.2)
# wmin = 0.45
# wmax = 1.2
# k = 0.3
# alpha = 3
# delta = 0.025
# phi_max = 0.09
# Lam = 12 * alpha * delta / (3 * alpha - 1)
# slope = -Lam * k / ((wmin + k * (Lam - 1)) * (wmax + k * (Lam - 1)))
# beta_0 = phi_max - slope * wmin
# theo = beta_0 + slope * w_range
# ax.plot(w_range, theo, '--')
for i, (g, d) in enumerate(pred_phi[pred_phi['interval'] != 'median'].groupby(['interval'], sort=False)):
    ax.fill_between(d['width'], d['lower'], d['upper'],
                    color=band_cor[g])
# ax.set_ylim([0.45, 1])
ax.set_xlim([0.5, 1.1])
ax.set_ylim([0, 0.1])

ax.set_ylabel('periplasmic biomass fraction\n$\phi_M$', fontsize=6)
ax.set_xlabel('w\naverage width [µm]', fontsize=6)
# plt.savefig('../../figures/Fig2_width_prediction_oe.pdf')
