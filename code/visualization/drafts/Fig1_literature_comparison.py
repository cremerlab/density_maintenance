# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
cor, pal = size.viz.matplotlib_style()
np.random.seed(666)
err_widths = {'95%': 0.25, '75%': 1, '25%': 2}
# %%
# Load Literature datasets
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
lit_mass_spec = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')
lit_size_data['aspect_ratio'] = lit_size_data['length_um'] / \
    lit_size_data['width_um']

lit_mass_spec = lit_mass_spec[lit_mass_spec['periplasm'] == True]
lit_mass_spec = lit_mass_spec.groupby(['dataset_name', 'condition',
                                      'growth_rate_hr']
                                      )['mass_frac'].sum().reset_index()
mass_spec_percs = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')
# %%
# Load our inferred data

# datasets
growth_rates = pd.read_csv(
    '../../data/mcmc/growth_parameter_percentiles.csv')
growth_rates = growth_rates[(growth_rates['strain'] == 'wildtype') &
                            (growth_rates['overexpression'] == 'none') &
                            (growth_rates['inducer_conc'] == 0)]

sizes = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
sizes = sizes[(sizes['strain'] == 'wildtype') &
              (sizes['overexpression'] == 'none') &
              (sizes['inducer_conc'] == 0)]
# sizes=sizes.groupby(['carbon_source'])[['width_median', 'length',
#                                           'volume', 'surface_to_volume', 'aspect_ratio']].agg(('mean', 'sem'))

# prot_data = pd.read_csv(
#     '../../data/summaries/summarized_protein_measurements.csv')
# prot_data = prot_data[(prot_data['strain'] == 'wildtype') &
#                       (prot_data['overexpression'] == 'none') &
#                       (prot_data['inducer_conc_ng_mL'] == 0)]
# prot_data = prot_data.groupby(['carbon_source'])[
#     'mass_frac'].agg(('mean', 'sem')).reset_index()

# Load parameters from mcmc
singulars = pd.read_csv('../../data/mcmc/growth_rate_linear_relations.csv')
# singulars = singulars[singulars['interval'] == 'median']

# Define markers and colors
markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h',
           'p', 'P', '*', 'o', '8', 'd', '>', 'v', '<', '^']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(lit_size_data['source'].unique())
for n in lit_mass_spec['dataset_name'].unique():
    names.append(n)
mapper = {n: {'m': m, 'c': c} for n, m, c in zip(names, markers, cors)}

# %%
# Make the plots of the widths lengths and volumes
fig, ax = plt.subplots(4, 2, figsize=(4, 6.6))
ax = ax.ravel()
for a in ax[:5]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[5].set_xlabel('surface area to volume [µm$^{-1}$]', fontsize=6)
ax[0].set_ylabel('average length [µm]', fontsize=6)
ax[1].set_ylabel('average width [µm]', fontsize=6)
ax[2].set_ylabel('average volume [µm$^{3}$]', fontsize=6)
ax[3].set_ylabel('surface area to volume [µm$^{-1}$]', fontsize=6)
ax[4].set_ylabel('average aspect ratio', fontsize=6)
ax[5].set_ylabel('periplasmic protein\nbiomass fraction', fontsize=6)
ax[-1].axis('off')
ax[-2].axis('off')
ax[0].set_ylim([1, 6])
ax[1].set_ylim([0.4, 1.3])
ax[2].set_ylim([0, 4.5])
ax[4].set_ylim([0, 8])
ax[5].set_ylim([0, 0.10])

# Plot the fits.
interval_colors = {'95%': cor['pale_blue'], '75%': cor['light_blue'],
                   '25%': cor['primary_blue'], 'median': cor['blue']}
lam_range = np.linspace(0, 2.5, 300)
axes = {'length': ax[0], 'width': ax[1], 'sav': ax[3],
        'alpha': ax[4], 'mass_fraction': ax[5]}
for i, (g, d) in enumerate(singulars.groupby(['quantity', 'interval'], sort=False)):
    if g[0] == 'mass_fraction':
        continue

    if g[1] == 'median':
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'],
                        '-', color=interval_colors[g[1]])
    else:
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'],
                                d['upper'], color=interval_colors[g[1]])


# Plot the literature data
alpha = 0.35
for g, d in lit_size_data.groupby(['source']):
    for i, v in enumerate(['length_um', 'width_um', 'volume_um3', 'surface_to_volume', 'aspect_ratio']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], linestyle='none', marker=mapper[g]['m'],
                   ms=3, markeredgecolor=cor['primary_black'],
                   markerfacecolor=mapper[g]['c'], alpha=alpha, label=g)
    ax[-1].plot([], [], linestyle='none', marker=mapper[g]['m'],
                ms=3, markeredgecolor=cor['primary_black'],
                markerfacecolor=mapper[g]['c'], alpha=alpha, label=g)

# for g, d in lit_mass_spec.groupby(['dataset_name']):
#     ax[5].plot(d['growth_rate_hr'], d['mass_frac'], linestyle='none', marker=mapper[g]['m'],
#                ms=3, markeredgecolor=cor['primary_black'],
#                markerfacecolor=mapper[g]['c'], alpha=alpha, label=g)
#     ax[-1].plot([], [], linestyle='none', marker=mapper[g]['m'],
#                 ms=3, markeredgecolor=cor['primary_black'],
#                 markerfacecolor=mapper[g]['c'], alpha=alpha, label=g)

for g, d in mass_spec_percs.groupby(['dataset_name']):
    meds = d[d['interval'] == 'median']
    ax[5].plot(meds[meds['quantity'] == 'mass_spec_sav']['lower'],
               meds[meds['quantity'] == 'mass_spec_phi_M']['lower'],
               marker=mapper[g]['m'], markeredgecolor=cor['primary_black'],
               markerfacecolor=mapper[g]['c'], alpha=alpha,
               linestyle='none', ms=3)

# Plot our data
for g, d in sizes.groupby(['carbon_source']):
    if (g == 'ezMOPS'):
        continue
    lam = growth_rates[(growth_rates['carbon_source'] == g)]
    lam_median = lam[lam['interval'] == 'median']
    for i, v in enumerate(['length_mu', 'width_mu', 'volume_mu', 'surface_area_vol_mu', 'aspect_ratio_mu']):
        _d = d[(d['quantity'] == v)]
        _d_median = _d[_d['interval'] == 'median']
        _d_cred = _d[_d['interval'] != 'median']
        # Add credible region
        for j, (__g, __d) in enumerate(lam[lam['interval'] != 'median'].groupby(['interval'], sort=False)):
            ax[i].hlines(_d_median['lower'], __d['lower'], __d['upper'],
                         color=cor['blue'], lw=err_widths[__g])
            ax[i].vlines(lam_median['lower'], _d_cred[_d_cred['interval'] == __g]['lower'], _d_cred[_d_cred['interval']
                         == __g]['upper'], color=cor['blue'], lw=err_widths[__g])

        ax[i].plot(lam_median['lower'], _d_median['lower'],
                   marker='o', markeredgecolor=cor['blue'],
                   markerfacecolor='white', markeredgewidth=0.75, ms=3, lw=1,
                   label='__nolegend__', color=cor['blue'])


for g, d in sizes.groupby(['carbon_source']):
    if (g == 'LB') | (g == 'ezMOPS'):
        continue
    med = d[d['interval'] == 'median']
    phiM_med = med[med['quantity'] == 'phi_M']['lower']
    sav_med = med[med['quantity'] == 'surface_area_vol_mu']['lower']
    for _g, _d in d[d['interval'] != 'median'].groupby(['interval'], sort=False):
        phi_M = _d[_d['quantity'] == 'phi_M']
        sav = _d[_d['quantity'] == 'surface_area_vol_mu']
        ax[5].vlines(sav_med, phi_M['lower'], phi_M['upper'], color=cor['blue'],
                     lw=err_widths[_g])
        ax[5].hlines(phiM_med, sav['lower'], sav['upper'], color=cor['blue'],
                     lw=err_widths[_g])
    ax[5].plot(sav_med, phiM_med, 'o', ms=3, markeredgewidth=0.75, markerfacecolor='white',
               markeredgecolor=cor['blue'], label='__nolegend__')

    # ax[5].plot(lam['lower'], [''], xerr=lam['sem'], yerr=d['sem'], lw=1,
    #    marker='o', linestyle='none', markeredgewidth=0.75, ms=4,
    #    markerfacecolor='white', markeredgecolor=cor['primary_blue'],
    #    label='__nolegend__', color=cor['blue'], capsize=0)

# Plot the parameter estimates
width = singulars[singulars['quantity'] == 'width_um']
ax[1].plot(width['growth_rate_hr'], width['lower'],
           '-', lw=1, color=cor['primary_blue'])
ax[0].plot(width['growth_rate_hr'], 3.3 * width['lower'],
           '--', lw=1, color=cor['primary_blue'])
ax[2].plot(width['growth_rate_hr'], (np.pi / 12) * width['lower'].values**2 * (3 * 3.3 *
           width['lower'].values - width['lower'].values), '--', lw=1, color=cor['primary_blue'])

# ax[-1].plot([], [], 'o', markeredgecolor=cor['primary_blue'],
#             markerfacecolor='white', markeredgewidth=0.5, ms=4, label='This study')
# ax[-1].legend()
# plt.tight_layout()
plt.savefig('../../figures/Fig1_literature_comparison.pdf',
            bbox_inches='tight')
# plt.savefig('../../figures/Fig1_literature_comparison_dataonly.pdf',
# bbox_inches='tight')
