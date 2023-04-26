# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
ppc_cmap = {'95%': cor['pale_black'], '75%': cor['light_black'],
            '50%': cor['primary_black'], '25%': cor['black'], '10%': cor['dark_black']}

# Load the various datasets
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
kdes = pd.read_csv('../../data/mcmc/literature_model_params_kde.csv')

size_data['aspect_ratio'] = size_data['length_um'].values / \
    size_data['width_um'].values
m2_ppcs = ppcs[ppcs['model'] == 'const_rho_mem']
ppcs = ppcs[ppcs['model'] == 'const_phi_mem']
kdes = kdes[kdes['model'] == 'const_phi_mem']
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']


# %%

# Instantiate figure canvas
fig, ax = plt.subplots(2, 2, figsize=(3.75, 2.5), sharex=True)
ax = ax.ravel()
ax[0].set_xlim([0, 2.25])
ax[0].set_ylim([0, 20])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 5])
ax[3].set_ylim([1, 8])

# Add titles and labels
ax[0].set_title('periplasmic protein mass per cell', fontsize=6)
ax[1].set_title('allocation towards membrane protein', fontsize=6)
ax[3].set_title('average cell aspect ratio', fontsize=6)
ax[2].set_title('membrane protein areal density', fontsize=6)

ax[0].set_ylabel('$m_{peri}$\n[fg / cell]', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$', fontsize=6)
ax[2].set_ylabel(r'$\rho_{mem}$'+'\n[fg / µm$^2$]', fontsize=6)
ax[3].set_ylabel(r'$\alpha$', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)


# Isolate the median values for the constants
const_medians = ppcs[(ppcs['interval'] == 'median') &
                     (ppcs['quantity'].isin(['alpha_sim', 'm_peri_sim',
                                             'phi_mem_sim']))]
const_rho_mem = m2_ppcs[(m2_ppcs['interval'] == 'median') &
                        (m2_ppcs['quantity'] == 'rho_mem_sim')]
# Plot the literature data
for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby(['dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['mass_fg'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=cor['light_purple'],
               markeredgewidth=0.5, alpha=0.75)

for g, d in ms_data[ms_data['localization'] == 'membrane'].groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=cor['pale_blue'],
               markeredgewidth=0.5, alpha=0.75)
    ax[2].plot(d['growth_rate_hr'], d['mass_fg'] / (2 * d['surface_to_volume'] * d['volume']),
               mapper[g]['m'], ms=3, markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5, alpha=0.75, color=cor['pale_blue'])

for g, d in size_data.groupby(['source']):
    ax[3].plot(d['growth_rate_hr'], d['aspect_ratio'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               markeredgewidth=0.5, alpha=0.75)

axes = {'m_peri_sim': ax[0], 'phi_mem_sim': ax[1], 'alpha_sim': ax[3]}
for i, (g, d) in enumerate(const_medians.groupby(['quantity'], sort=False)):
    d.sort_values(by='growth_rate_hr')
    axes[g].plot(d['growth_rate_hr'], d['lower'], '-',
                 lw=0.75, color=cor['primary_black'])

ax[2].plot(const_rho_mem['growth_rate_hr'], const_rho_mem['lower'], '-',
           lw=0.75, color=cor['primary_black'])
plt.subplots_adjust(hspace=0.35)
plt.savefig('../../Fig3_constant_plots.pdf', bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 3, figsize=(3, 1))
ax = ax.ravel()
axes = {'m_peri_mu': ax[0], 'phi_mem_mu': ax[1],
        'alpha': ax[2], 'rho_mem_mu': ax[3]}
for a in ax:
    a.set_yticks([])
    a.set_facecolor('#FFF')
ax[0].set_xlim([8, 12])
ax[0].set_xticks([8, 10, 12])
ax[1].set_xlim([0.10, 0.15])
ax[1].set_xticks([0.10, 0.125, 0.15])
ax[2].set_xlim([3.5, 4.5])
ax[0].set_xlabel('$m_{peri}$ [fg / cell]', fontsize=6)
ax[1].set_xlabel('$\phi_{mem}$', fontsize=6)
ax[3].set_ylim([0, 5])

# Plot the probability density
percs = [(2.5, 97.5), (12.5, 87.5), (37.5, 62.5), (49.5, 50.5)]
perc_colors = {'95%': 'pale',
               '75%': 'light',
               '25%': 'primary',
               'median': 'dark'}
for g, d in kdes[kdes['parameter'].isin(['m_peri_mu', 'phi_mem_mu', 'alpha', 'rho_mem_mu'])].groupby(['parameter']):
    axes[g].fill_between(d['value'], np.zeros(len(d)), d['kde'] / d['kde'].max(),  color=cor['pale_blue'],
                         alpha=0.5)
    axes[g].plot(d['value'], d['kde'] / d['kde'].max(),
                 '-', lw=0.5, color=cor['primary_blue'])

for i, (g, d) in enumerate(m2_ppcs[m2_ppcs['quantity'].isin(['w_rep', 'ell_rep', 'vol_rep'])].groupby(['quantity', 'interval'], sort=False)):
    if g[1] != 'median':
        ax[i].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                           color=cor[f'{perc_colors[g[0]]}_green'], alpha=0.5)
# # Plot the major percentiles
# axes = {'m_peri_sim': ax[0], 'phi_mem_sim': ax[1], 'alpha_sim': ax[2]}
# err_widths = {'95%': 0.5, '75%': 1, '25%': 1.5}
# for g, d in ppcs[ppcs['quantity'].isin(['alpha_sim', 'phi_mem_sim', 'alpha_sim'])].groupby(['quantity', 'interval']):
#     if g[1] in ['95%', '75%', '25%', 'median']:
#         if g[1] != 'median':
#             axes[g[0]].hlines(-0.05, d['lower'].values[0], d['upper'].values[0],
#                               lw=err_widths[g[1]], color=cor['primary_blue'])

# plt.tight_layout()
plt.savefig('../../figures/Fig3_const_param_posteriors.pdf',
            bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 4, figsize=(6, 1.5), sharex=True)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_xlim([0, 2.6])
ax[1].set_ylabel('width [µm]', fontsize=6)
ax[2].set_ylabel('length [µm]', fontsize=6)
ax[3].set_ylabel('volume [µm$^3$]', fontsize=6)
ax[0].set_ylabel('aspect ratio', fontsize=6)
ax[1].set_ylim([0.45, 1.25])
ax[2].set_ylim([1, 6.5])
ax[3].set_ylim([-0.05, 5])
ax[0].set_ylim([1, 8])
perc_colors = {'95%': 'pale',
               '75%': 'light',
               '25%': 'primary',
               'median': 'dark'}
for g, d in size_data.groupby(['source']):
    for i, p in enumerate(['aspect_ratio', 'width_um', 'length_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[p], mapper[g]['m'], ms=4,
                   markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
                   alpha=0.75, zorder=1000)
axes = {'w_rep': ax[1], 'ell_rep': ax[2], 'vol_rep': ax[3], 'alpha_rep': ax[0]}
for g, d in ppcs[ppcs['quantity'].isin(['w_rep', 'ell_rep', 'vol_rep', 'alpha_rep']) &
                 (ppcs['interval'].isin(perc_colors.keys()))].groupby(['quantity', 'interval'], sort=False):
    if g[1] != 'median':
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
    else:
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'], lw=1,
                        color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
plt.tight_layout()
plt.savefig('../../figures/Fig3_size_ppc.pdf')
# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 1.5), sharex=True)
for a in ax:
    a.set_xlabel('average cell width [µm]', fontsize=6)
    a.set_ylim([0, 1])
    a.set_xlim([0.55, 1.05])
    a.set_ylabel(r'$\phi_{peri} / \phi_{mem}$' +
                 '\nrelative envelope allocation', fontsize=6)
for g, d in ms_data[ms_data['localization'].isin(['periplasm', 'membrane'])].groupby(['dataset_name']):
    peri = d[d['localization'] == 'periplasm']
    mem = d[d['localization'] == 'membrane']
    for a in ax:
        a.plot(peri['width'], peri['mass_frac'].values/mem['mass_frac'].values, mapper[g]['m'],
               color=mapper[g]['c'], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
               alpha=0.5, ms=5)

for g, d in ppcs[ppcs['interval'].isin(perc_colors.keys()) & ppcs['quantity'].isin(['rel_phi_rep'])].groupby(['interval'], sort=False):
    if g[0] != 'median':
        ax[0].fill_between(d['width'], d['lower'], d['upper'],
                           alpha=0.3, color=cor[f'{perc_colors[g]}_blue'])
    else:
        ax[0].plot(d['width'], d['lower'], d['upper'],
                   alpha=0.3, color=cor[f'{perc_colors[g]}_blue'])

for g, d in m2_ppcs[m2_ppcs['interval'].isin(perc_colors.keys()) & m2_ppcs['quantity'].isin(['rel_phi_rep'])].groupby(['interval'], sort=False):
    if g[0] != 'median':
        ax[1].fill_between(d['width'], d['lower'], d['upper'],
                           alpha=0.3, color=cor[f'{perc_colors[g]}_green'])
    else:
        ax[1].plot(d['width'], d['lower'], d['upper'],
                   alpha=0.3, color=cor[f'{perc_colors[g]}_green'])

plt.savefig('../../figures/Fig3_model_comparison.pdf')
# %%
fig, ax = plt.subplots(2, 1, figsize=(2, 2.5), sharex=True)
ax[0].set_ylabel(r'$\rho_{mem}$ [fg / µm$^2$]', fontsize=6)
ax[1].set_ylabel(r'$\rho_{peri}$ [fg / µm$^3$]', fontsize=6)
ax[0].set_ylim([0, 5])
ax[1].set_ylim([0, 150])
ax[1].set_xlim([0, 2])
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
axes = {'membrane': ax[0], 'periplasm': ax[1],
        'rho_mem_rep': ax[0], 'rho_peri_rep': ax[1]}
ms_colors = {'membrane': cor['pale_blue'], 'periplasm': cor['light_purple']}

for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm'])].groupby(['dataset_name', 'localization']):
    if g[1] == 'membrane':
        d['density'] = d['mass_fg'] / \
            (2 * d['surface_to_volume'] * d['volume'])
    else:
        d['density'] = d['mass_fg'] / \
            (d['surface_to_volume'] * d['volume'] * 0.0249)

    axes[g[1]].plot(d['growth_rate_hr'], d['density'], mapper[g[0]]['m'], markeredgecolor=cor['primary_black'],
                    color=ms_colors[g[1]], alpha=0.75, ms=4)

perc_colors = {'rho_peri_rep': {'95%': cor['pale_purple'], '75%': cor['light_purple'], '25%': cor['primary_purple'], 'median': cor['purple']},
               'rho_mem_rep': {'95%': cor['pale_blue'], '75%': cor['light_blue'], '25%': cor['primary_blue'], 'median': cor['blue']}}
for g, d in ppcs[ppcs['quantity'].isin(['rho_mem_rep', 'rho_peri_rep'])].groupby(['quantity', 'interval']):
    if g[1] in ['95%', '75%', '25%', 'median']:
        if g[1] != 'median':
            axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                    color=perc_colors[g[0]][g[1]], alpha=0.5)
        else:
            axes[g[0]].plot(d['growth_rate_hr'], d['lower'], '-',
                            lw=1, color=perc_colors[g[0]][g[1]], alpha=0.5)
plt.savefig('../../figures/Fig3_density_ppc.pdf')

# %%
# Compile the relative allocation from the ms data
rel_ms = ms_data[ms_data['localization'].isin(['periplasm', 'membrane'])]

fig, ax = plt.subplots(1, 1, figsize=(3, 2.5))
ax.set_ylim([0, 1])
ax.set_xlim([0.55, 1])
ax.set_ylabel(
    '$\phi_{peri} / \phi_{mem}$\nrelative allocation in envelope', fontsize=6)
ax.set_xlabel('average cell width [µm]', fontsize=6)
perc_colors = {'95%': cor['pale_black'], '75%': cor['light_black'],
               '25%': cor['primary_black'], 'median': cor['black']}

for g, d in ppcs[ppcs['quantity'] == 'rel_phi_rep'].groupby(['interval'], sort=False):
    if g in perc_colors.keys():
        if g != 'median':
            ax.fill_between(d['width'].values, d['lower'].values,
                            d['upper'].values, color=perc_colors[g], alpha=0.5)
        else:
            ax.plot(d['width'].values, d['lower'].values,
                    '-', lw=1, color=perc_colors[g])

for g, d in rel_ms.groupby(['dataset_name']):
    mem = d[d['localization'] == 'membrane']
    peri = d[d['localization'] == 'periplasm']
    ax.plot(mem['width'],  peri['mass_frac'].values / mem['mass_frac'].values, mapper[g]['m'],
            color=mapper[g]['c'], ms=4, markeredgecolor=cor['primary_black'],
            markeredgewidth=0.5, alpha=0.75)
plt.savefig('../../figures/Fig3_relative_allocation_prediction.pdf')
