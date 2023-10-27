# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

prot = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
dna = pd.read_csv('../../data/literature/collated_dna_protein_ratio.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
drymass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')
rho_cyt = pd.read_csv(
    '../../data/mcmc/approximated_cytoplasmic_density_wide.csv')
ppcs = pd.read_csv('../../data/mcmc/approximated_cytoplasmic_density_ppcs.csv')
fig, ax = plt.subplots(1, 1, figsize=(3, 3))
ax.set_ylim(100, 700)
for g, d in drymass.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)

for g, d in rho_cyt.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['primary_green']
    ax.vlines(d['growth_rate_hr'], d['2.5%'],
              d['97.5%'], lw=0.5, color=cor['primary_green'])
    ax.plot(d['growth_rate_hr'], d['median_value'], **fmt, zorder=1000)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('mass density [fg / fL]', fontsize=6)
# ax.legend()
plt.savefig('../../figures/FigA1_cytoplasmic_drymass_density.pdf',
            bbox_inches='tight')

# %%
fig, ax = plt.subplots(2, 2, figsize=(3, 3), sharex=True)

for i in range(2):
    ax[1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax = ax.ravel()
for a in ax:
    a.set_xlim([0, 2.5])
ax[0].set_ylim([100, 800])
ax[1].set_ylim(0, 0.1)
ax[2].set_ylim([0, 3])
ax[3].set_ylim([0, 10])
ax[0].set_ylabel('$M_{prot}^{(tot)}$ [fg / cell]\ntotal protein', fontsize=6)
ax[1].set_ylabel(r'$\theta_{DNA}$' + '\nDNA-to-protein', fontsize=6)
ax[2].set_ylabel('V [fL]\naverage volume', fontsize=6)
ax[3].set_ylabel('$S_A$ [µm$^2$]\naverage surface area', fontsize=6)

for g, d in prot.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['primary_red']
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt, ms=4)

for g, d in dna.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['primary_blue']
    ax[1].plot(d['growth_rate_hr'], d['DNA_protein_ratio'], **fmt, ms=4)

for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['primary_black']
    ax[2].plot(d['growth_rate_hr'], d['volume_um3'], **fmt, ms=4)
    fmt['color'] = cor['primary_purple']
    ax[3].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt, ms=4)


axes = {'M_prot_pred': [ax[0], 'red'], 'DNA_pred': [ax[1], 'blue'],
        'vol_pred': [ax[2], 'black'], 'SA_pred': [ax[3], 'purple']}
inter_colors = {'95%': 'pale_', '75%': 'light_',
                '25%': 'primary_', 'median': ''}

for g, d in ppcs.groupby('quantity'):
    a = axes[g][0]
    for _g, _d in d.groupby('interval', sort=False):
        a.fill_between(_d['growth_rate_hr'], _d['lower'], _d['upper'],
                       color=cor[f'{inter_colors[_g]}{axes[g][1]}'], alpha=0.75)

for a in ax:
    a.legend(fontsize=5)
plt.savefig('../../figures/FigA1_empirical_fits.pdf')
# %%
