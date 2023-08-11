# %%
import scipy.stats
import numpy as np
import pandas as pd
import cmdstanpy
import size.viz
import matplotlib.pyplot as plt
import arviz as az
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./width_scaling_periplasm.stan')


# %%
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
biomass_data = pd.read_csv(
    '../../data/literature/collated_drymass_densities.csv')
prot_per_cell = pd.read_csv(
    '../../data/literature/collated_protein_per_cell.csv')
dna_ratio = pd.read_csv('../../data/literature/collated_dna_protein_ratio.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']

mem_ms = ms_data[ms_data['localization'] == 'membrane'].copy()
mem_ms['rho'] = mem_ms['mass_fg'] / (2 * mem_ms['surface_area'])

width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])

size_data['phiRb'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width'] = width_popt[1] + \
    width_popt[0] * phiRb_data['growth_rate_hr']


# %%
N_ppc = 200
lam_range = np.linspace(0, 2.5, N_ppc)
phiRb_range = phiRb_popt[1] + phiRb_popt[0] * lam_range
data_dict = {
    'N_prot': len(prot_per_cell),
    'N_phiRb': len(phiRb_data),
    'N_ms': len(mem_ms),
    'N_biomass': len(biomass_data),
    'N_size': len(size_data),
    'N_dna': len(dna_ratio),
    'N_ppc': N_ppc,

    'prot_per_cell': prot_per_cell['fg_protein_per_cell'].values,
    'prot_per_cell_lam': prot_per_cell['growth_rate_hr'].values,

    'dna_to_protein': dna_ratio['DNA_protein_ratio'].values,
    'dna_lam': dna_ratio['growth_rate_hr'].values,
    'phiRb': phiRb_data['mass_fraction'].values,
    'phiRb_lam': phiRb_data['growth_rate_hr'].values,

    'phi_mem': mem_ms['mass_frac'].values,
    'rho_mem': mem_ms['rho'].values,
    'm_peri': ms_data[ms_data['localization'] == 'periplasm']['mass_fg'].values,

    'rho_biomass': biomass_data['drymass_density_fg_fL'].values,

    'aspect_ratio': size_data['length_um'].values / size_data['width_um'].values,
    'width': size_data['width_um'].values,
    'vol': size_data['volume_um3'].values,
    'size_lam': size_data['growth_rate_hr'].values,

    'lam_range': lam_range,
    'phiRb_range': phiRb_range

}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
ax[0].set_xlabel('ribosomal allocation $\phi_{Rb}$', fontsize=6)
ax[0].set_ylabel('average width [µm]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_ylabel('average volume [µm$^3$]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[2].set_ylabel(r'periplasm allocation $\phi_{peri}$', fontsize=6)

uppers = [97.5, 87.5, 62.5, 50]
lowers = [2.5, 12.5, 37.5, 50]
labels = ['95%', '75%', '25%', 'median']


pars = ['width_pred_sim', 'vol_pred_sim', 'phi_peri_sim']
reps = []
for i, p in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    for j in range(len(lam_range)):
        post.loc[post[f'{p}_dim_0'] == j, 'lam'] = lam_range[j]
        post.loc[post[f'{p}_dim_0'] == j, 'phiRb'] = phiRb_range[j]
    percs = size.viz.compute_percentiles(post, p, ['lam', 'phiRb'],
                                         lower_bounds=lowers,
                                         upper_bounds=uppers,
                                         interval_labels=labels)
    reps.append(percs)


uppers = [97.5, 87.5, 62.5, 50]
lowers = [2.5, 12.5, 37.5, 50]
labels = ['95%', '75%', '25%', 'median']
perc_colors = {'95%': cor['pale_blue'], '75%': cor['light_blue'],
               '25%': cor['primary_blue'], 'median': cor['blue']}
for i, r in enumerate(reps):
    if i == 0:
        x = 'phiRb'
    else:
        x = 'lam'
    for g, d in r.groupby('interval', sort=False):
        ax[i].fill_between(d[x], d['lower'], d['upper'],
                           alpha=0.5, color=perc_colors[g])

for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['phiRb'], d['width_um'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)

for g, d in phiRb_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['mass_fraction'], d['width'], **fmt)

for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], **fmt)

plt.savefig('./scaling_with_periplasm.pdf')

# %%
fig, ax = plt.subplots(2, 3, figsize=(7, 4))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylim([100, 400])
ax[0].set_ylabel(r'biomass density $\rho_{bio}$ [fg / µm$^3$]', fontsize=6)
ax[1].set_ylim([0, 6])
ax[1].set_ylabel(r'membrane density $\rho_{mem}$ [fg / µm$^2$]', fontsize=6)
ax[2].set_ylim([0, 0.2])
ax[2].set_ylabel(r'membrane allocation $\phi_{mem}$', fontsize=6)
ax[3].set_ylim([0, 30])
ax[3].set_ylabel('periplasm mass $m_{peri}$ [fg/cell]', fontsize=6)
ax[3].set_xlim([0, 1.5])
ax[4].set_ylim([0, 0.1])
ax[4].set_ylabel('DNA to total protein ratio', fontsize=6)
ax[5].axis('off')
uppers = [97.5, 87.5, 62.5, 50]
lowers = [2.5, 12.5, 37.5, 50]
labels = ['95%', '75%', '25%', 'median']
pars = ['rho_biomass_rep', 'rho_mem_rep',
        'phi_mem_rep', 'm_peri_rep', 'dna_to_protein_rep']
perc_colors = {'95%': cor['pale_black'], '75%': cor['light_black'],
               '25%': cor['primary_black'], 'median': cor['dark_black']}
for i, p in enumerate(pars):
    _df = samples.posterior[p].to_dataframe().reset_index()
    _df['idx'] = 1
    _percs = size.viz.compute_percentiles(_df, p, 'idx',
                                          lower_bounds=lowers,
                                          upper_bounds=uppers,
                                          interval_labels=labels)
    for g, d in _percs.groupby('interval'):
        if g == 'median':
            ax[i].plot(lam_range, d['lower'].values[0] * np.ones_like(lam_range),  alpha=0.5,
                       lw=1, color=perc_colors[g])
        else:
            ax[i].fill_between(lam_range, d['lower'], d['upper'],
                               alpha=0.5, color=perc_colors[g])

for g, d in biomass_data.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.75)
    ax[0].plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)
    ax[5].plot([], [], **fmt)

for g, d in mem_ms.groupby('dataset_name'):
    fmt = size.viz.style_point(g, alpha=0.75)
    ax[1].plot(d['growth_rate_hr'], d['rho'], **fmt)
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], **fmt)
    ax[5].plot([], [], **fmt)

for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby('dataset_name'):
    fmt = size.viz.style_point(g, alpha=0.75)
    ax[3].plot(d['growth_rate_hr'], d['mass_fg'], **fmt)
    ax[5].plot([], [], **fmt)

for g, d in dna_ratio.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.75)
    ax[4].plot(d['growth_rate_hr'], d['DNA_protein_ratio'], **fmt)
    ax[5].plot([], [], **fmt)

ax[5].legend(ncols=2, bbox_to_anchor=(1.5, 1), fontsize=5, handletextpad=0.5)
plt.tight_layout()
plt.savefig('./constants.pdf')
