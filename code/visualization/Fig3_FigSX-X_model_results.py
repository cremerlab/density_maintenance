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
ppcs = pd.read_csv(
    '../../data/mcmc/literature_model_params.csv')

m1_ppcs = ppcs[(ppcs['volume_scale'] == 'linear_width')
               & (ppcs['model'] == 'const_phi_mem')]
m1_ppcs_smk = ppcs[(ppcs['volume_scale'] == 'smk')
                   & (ppcs['model'] == 'const_phi_mem')]
m2_ppcs = ppcs[(ppcs['volume_scale'] == 'linear_width')
               & (ppcs['model'] == 'const_rho_mem')]
size_data['aspect_ratio'] = size_data['length_um'].values / \
    size_data['width_um'].values
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']

# %%

# Instantiate figure canvas
fig, ax = plt.subplots(2, 2, figsize=(3.95, 2.5), sharex=True)
ax = ax.ravel()
ax[3].set_xlim([0, 2.25])
ax[3].set_ylim([0, 20])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 5])
ax[0].set_ylim([1, 8])

# Add titles and labels
ax[0].set_title('(i) average cell aspect ratio', fontsize=6)
ax[1].set_title('(ii) allocation towards membrane protein', fontsize=6)
ax[2].set_title('(iii) membrane protein areal density', fontsize=6)
ax[3].set_title('(iv) periplasmic protein mass per cell', fontsize=6)

ax[3].set_ylabel('$m_{peri}$\n[fg / cell]', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$', fontsize=6)
ax[2].set_ylabel(r'$\rho_{mem}$'+'\n[fg / µm$^2$]', fontsize=6)
ax[0].set_ylabel(r'$\alpha$', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)


# Isolate the median values for the constants
const_medians = m1_ppcs[(m1_ppcs['interval'] == 'median') &
                        (m1_ppcs['quantity'].isin(['alpha_sim', 'm_peri_sim',
                                                   'phi_mem_sim']))]
const_rho_mem = m2_ppcs[(m2_ppcs['interval'] == 'median') &
                        (m2_ppcs['quantity'] == 'rho_mem_sim')]
# Plot the literature data
for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby(['dataset_name']):
    ax[3].plot(d['growth_rate_hr'], d['mass_fg'], mapper[g]['m'],
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
    ax[0].plot(d['growth_rate_hr'], d['aspect_ratio'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               markeredgewidth=0.5, alpha=0.75)

axes = {'m_peri_sim': ax[3], 'phi_mem_sim': ax[1], 'alpha_sim': ax[0]}
for i, (g, d) in enumerate(const_medians.groupby(['quantity'], sort=False)):
    d.sort_values(by='growth_rate_hr')
    axes[g].plot(d['growth_rate_hr'], d['lower'], '-',
                 lw=0.75, color=cor['primary_black'])

ax[2].plot(const_rho_mem['growth_rate_hr'], const_rho_mem['lower'], '-',
           lw=0.75, color=cor['primary_black'])
plt.subplots_adjust(hspace=0.35)
plt.savefig('../../Fig3_constant_plots.pdf', bbox_inches='tight')

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


for g, d in m1_ppcs[m1_ppcs['quantity'].isin(['w_rep', 'ell_rep', 'vol_rep', 'alpha_rep']) &
                    (m1_ppcs['interval'].isin(perc_colors.keys()))].groupby(['quantity', 'interval'], sort=False):
    if g[1] != 'median':
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
    else:
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'], lw=1,
                        color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
plt.tight_layout()
plt.savefig('../../figures/Fig3_size_ppc.pdf')

# %%
fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$m_{peri}$\n[fg / cell]', fontsize=6)
ax[0].set_title('periplasmic protein per cell', fontsize=6)
ax[0].set_xlim([0, 2.1])
ax[0].set_ylim([0, 20])
ax[1].set_title('(i) prediction from model I', fontsize=6)
ax[2].set_title('(ii) prediction from model II', fontsize=6)
for a in ax[1:]:
    a.set_xlabel('average cell width [µm]', fontsize=6)
    a.set_ylim([0, 1])
    a.set_xlim([0.55, 1.05])
    a.set_ylabel(r'$\phi_{peri} / \phi_{mem}$' +
                 '\nrelative envelope allocation', fontsize=6)

for g, d in ms_data[ms_data['localization'].isin(['periplasm', 'membrane'])].groupby(['dataset_name']):
    peri = d[d['localization'] == 'periplasm']
    mem = d[d['localization'] == 'membrane']
    for a in ax[1:]:
        a.plot(peri['width'], peri['mass_frac'].values/mem['mass_frac'].values, mapper[g]['m'],
               color=mapper[g]['c'], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
               alpha=0.5, ms=5)

for g, d in ms_data[ms_data['localization'].isin(['periplasm'])].groupby('dataset_name'):
    ax[0].plot(d['growth_rate_hr'], d['mass_fg'], mapper[g]['m'], ms=4,
               markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               alpha=0.74, zorder=1000)

for g, d in m1_ppcs[m1_ppcs['interval'].isin(perc_colors.keys()) & m1_ppcs['quantity'].isin(['m_peri_rep'])].groupby(['interval'], sort=False):
    if g[0] != 'median':
        ax[0].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                           alpha=0.3, color=cor[f'{perc_colors[g]}_black'])
    else:
        ax[0].plot(d['growth_rate_hr'], d['lower'], d['upper'],
                   alpha=0.3, color=cor[f'{perc_colors[g]}_black'])


for g, d in m1_ppcs[m1_ppcs['interval'].isin(perc_colors.keys()) & m1_ppcs['quantity'].isin(['rel_phi_rep'])].groupby(['interval'], sort=False):
    if g[0] != 'median':
        ax[1].fill_between(d['width'], d['lower'], d['upper'],
                           alpha=0.3, color=cor[f'{perc_colors[g]}_blue'])
    else:
        ax[1].plot(d['width'], d['lower'], d['upper'],
                   alpha=0.3, color=cor[f'{perc_colors[g]}_blue'])


for g, d in m2_ppcs[m2_ppcs['interval'].isin(perc_colors.keys()) & m2_ppcs['quantity'].isin(['rel_phi_rep'])].groupby(['interval'], sort=False):
    if g[0] != 'median':
        ax[2].fill_between(d['width'], d['lower'], d['upper'],
                           alpha=0.3, color=cor[f'{perc_colors[g]}_green'])
    else:
        ax[2].plot(d['width'], d['lower'], d['upper'],
                   alpha=0.3, color=cor[f'{perc_colors[g]}_green'])
plt.tight_layout()
plt.savefig('../../figures/Fig3_model_comparison.pdf')

# ##############################################################################
# ASSOCIATED SUPPLEMENTAL FIGURES
# ##############################################################################
# %%
fig, ax = plt.subplots(1, 4, figsize=(6, 1.5), sharex=True)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_xlim([0, 2.6])
ax[0].set_ylabel('aspect ratio', fontsize=6)
ax[1].set_ylabel('volume [µm$^3$]', fontsize=6)
ax[2].set_ylabel('width [µm]', fontsize=6)
ax[3].set_ylabel('length [µm]', fontsize=6)
ax[0].set_ylim([1, 8])
ax[1].set_ylim([-0.05, 5])
ax[2].set_ylim([0.45, 1.25])
ax[3].set_ylim([1, 6.5])
perc_colors = {'95%': 'pale',
               '75%': 'light',
               '25%': 'primary',
               'median': 'dark'}
for g, d in size_data.groupby(['source']):
    for i, p in enumerate(['aspect_ratio', 'volume_um3', 'width_um', 'length_um']):
        ax[i].plot(d['growth_rate_hr'], d[p], mapper[g]['m'], ms=4,
                   markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
                   alpha=0.75, zorder=1000)
axes = {'w_rep': ax[2], 'ell_rep': ax[3], 'vol_rep': ax[1], 'alpha_rep': ax[0]}


for g, d in m1_ppcs_smk[m1_ppcs_smk['quantity'].isin(['w_rep', 'ell_rep', 'vol_rep', 'alpha_rep']) &
                        (m1_ppcs_smk['interval'].isin(perc_colors.keys()))].groupby(['quantity', 'interval'], sort=False):
    if g[1] != 'median':
        axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
    else:
        axes[g[0]].plot(d['growth_rate_hr'], d['lower'], lw=1,
                        color=cor[f'{perc_colors[g[1]]}_black'], alpha=0.25)
plt.tight_layout()
plt.savefig('../../figures/FigSX_SMK_size_ppc.pdf')

# %%
fig, ax = plt.subplots(2, 3, figsize=(6, 3), sharex=True)

# Format and label axes
for i in range(3):
    ax[1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
for i in range(2):
    ax[i, 0].set_ylabel('$\phi_{mem}$', fontsize=6)
    ax[i, 1].set_ylabel(r'$\rho_{mem}$ [fg / $\mu$m$^2$]', fontsize=6)
    ax[i, 2].set_ylabel(r'$\rho_{peri}$ [fg / $\mu$m$^3$]', fontsize=6)
ax[0, 0].set_title('allocation towards membrane proteins', fontsize=6)
ax[0, 1].set_title('membrane protein areal density', fontsize=6)
ax[0, 2].set_title('periplasmic protein density', fontsize=6)

# Change limits
for a in ax.ravel():
    a.set_xlim([0, 2.1])
for i in range(2):
    ax[i, 0].set_ylim([0, 1])
    ax[i, 1].set_ylim([0, 5])
    ax[i, 2].set_ylim([0, 150])
# Add data plots
for g, d in ms_data.groupby(['dataset_name']):
    mem = d[d['localization'] == 'membrane']
    peri = d[d['localization'] == 'periplasm']
    for i in range(2):
        ax[i, 0].plot(mem['growth_rate_hr'], mem['mass_frac'], mapper[g]['m'],
                      markeredgecolor=cor['primary_black'], markeredgewidth=0.25,
                      color=mapper[g]['c'], ms=4, alpha=0.25, zorder=1000)
        ax[i, 1].plot(mem['growth_rate_hr'], mem['mass_fg'] / (2 * mem['surface_area']), mapper[g]['m'],
                      markeredgecolor=cor['primary_black'], markeredgewidth=0.25,
                      color=mapper[g]['c'], ms=4, alpha=0.25, zorder=1000)
        ax[i, 2].plot(peri['growth_rate_hr'], peri['mass_fg'] / (peri['surface_area'] * 0.025), mapper[g]['m'],
                      markeredgecolor=cor['primary_black'], markeredgewidth=0.25,
                      color=mapper[g]['c'], ms=4, alpha=0.25, zorder=1000)

# Add prediction bands
cols = {'phi_mem_rep': 0, 'rho_mem_rep': 1, 'rho_peri_rep': 2}
shades = ['blue', 'green']
for i, m in enumerate([m1_ppcs, m2_ppcs]):
    for g, d in m[m['quantity'].isin(['phi_mem_rep', 'rho_mem_rep', 'rho_peri_rep']) &
                  m['interval'].isin(perc_colors.keys())].groupby(['quantity', 'interval'], sort=False):

        ax[i, cols[g[0]]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                       color=cor[f'{perc_colors[g[1]]}_{shades[i]}'], alpha=0.75)
plt.tight_layout()
plt.savefig('../../figures/FigSX_density_model_comparison_ppc.pdf')
