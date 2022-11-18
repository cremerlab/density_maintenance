# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
cor, pal = size.viz.matplotlib_style()

# Load the various datasets
grdata = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')
mass_spec = pd.read_csv(
    '../../data/calculated_envelope_protein_mass_and_concentration.csv')
SAV_data = pd.read_csv('../../data/si2017_SAV.csv')
mass_spec['phi_peri'] = mass_spec['periplasmic_protein_mass_fg'].values / \
    mass_spec['total_protein_mass_fg']
sizedata = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
protdata = pd.read_csv(
    '../../data/protein_quantification/protein_quantification_summary.csv')
cmap = sns.color_palette("ch:start=.2,rot=-.3", n_colors=8).as_hex()
carbons = ['LB', 'ezMOPS', 'glucoseCAA',
           'glucose', 'glycerol', 'sorbitol', 'acetate']
carb_cor = {c: cmap[-i - 1] for i, c in enumerate(carbons)}

# %%
fig, ax = plt.subplots(2, 1, figsize=(2.5, 2.25), sharex=True)

# Plot the mass spec data
i = 0
marks = ['o', 's', 'v', '^']
for g, d in mass_spec.groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['phi_peri'], marks[i], color=cor['light_black'],
               alpha=0.5, ms=4, label=g)
    i += 1

ax[0].plot(SAV_data['growth_rate_hr'], SAV_data['SAV_inv_um'],
           'X', alpha=0.5, color=cor['light_black'], ms=4, label='Si et al. 2017')

for g, d in sizedata.groupby(['carbon_source']):
    if g not in carbons:
        continue
    _gr = grdata[grdata['carbon_source'] == g]
    _SAV = d[d['parameter'] == 'SAV_inv_um']
    _prot = protdata[(protdata['carbon_source'] == g) &
                     (protdata['parameter'] == 'phi_peri')]
    print(g)
    if (len(_gr) > 0) & (len(_prot) > 0):
        _color = carb_cor[g]

        ax[0].vlines(_gr['median'], _SAV['2.5%'],
                     _SAV['97.5%'], lw=0.5, color=_color)
        ax[0].vlines(_gr['median'], _SAV['12.5%'],
                     _SAV['87.5%'], lw=1.5, color=_color)
        ax[0].hlines(_SAV['median'], _gr['2.5%'],
                     _gr['97.5%'], lw=0.5, color=_color)
        ax[0].hlines(_SAV['median'], _gr['12.5%'],
                     _gr['87.5%'], lw=1.5, color=_color)
        ax[0].plot(_gr['median'], _SAV['median'],
                   'o', ms=3, markeredgewidth=0.5, markerfacecolor='w',
                   markeredgecolor=_color)
        if g != 'LB':
            # Protein data
            ax[1].vlines(_gr['median'], _prot['2.5%'],
                         _prot['97.5%'], lw=0.5, color=_color)
            ax[1].vlines(_gr['median'], _prot['12.5%'],
                         _prot['87.5%'], lw=1.5, color=_color)
            ax[1].hlines(_prot['median'], _gr['2.5%'],
                         _gr['97.5%'], lw=0.5, color=_color)
            ax[1].hlines(_prot['median'], _gr['12.5%'],
                         _gr['87.5%'], lw=1.5, color=_color)
            ax[1].plot(_gr['median'], _prot['median'],
                       'o', ms=3, markeredgewidth=0.5, markerfacecolor='w',
                       markeredgecolor=_color)


# Format axes
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel(r'$\frac{A_{surface}}{V_{cell}}$ [µm$^{-1}$]', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$', fontsize=6)

ax[1].set_ylim([0, 0.15])
ax[1].set_yticks([0, 0.05, 0.1, 0.15])
ax[1].legend(fontsize=4)
ax[0].legend(fontsize=4)
plt.tight_layout()
plt.savefig('../../figures/Fig1_SAV_mass_frac.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), sharex=True)
for g, d in sizedata.groupby(['carbon_source']):
    if g not in carbons:
        continue
    # _gr = grdata[grdata['carbon_source'] == g]
    _SAV = d[d['parameter'] == 'SAV_inv_um']
    _prot = protdata[(protdata['carbon_source'] == g) &
                     (protdata['parameter'] == 'phi_peri')]
    print(g)
    if (len(_gr) > 0) & (len(_prot) > 0):
        _color = carb_cor[g]
        ax.vlines(_SAV['median'], _prot['2.5%'],
                  _prot['97.5%'], lw=0.5, color=_color)
        ax.vlines(_SAV['median'], _prot['12.5%'],
                  _prot['87.5%'], lw=1.5, color=_color)
        ax.hlines(_prot['median'], _SAV['2.5%'],
                  _SAV['97.5%'], lw=0.5, color=_color)
        ax.hlines(_prot['median'], _SAV['12.5%'],
                  _SAV['87.5%'], lw=1.5, color=_color)
        ax.plot(_SAV['median'], _prot['median'],
                'o', ms=4, markeredgewidth=0.5, markerfacecolor='w',
                markeredgecolor=_color)


ax.set_xlim([5.5, 7.75])
ax.set_ylim([0, 0.15])
ax.set_xticks([6, 6.5, 7, 7.5])
ax.set_yticks([0, 0.05, 0.1, 0.15])
ax.set_xlabel('surface area to volume [µm$^{-1}$]\n$Q_{AV}$', fontsize=6)
ax.set_ylabel('$\phi_{peri}$\nperiplasmic mass fraction', fontsize=6)
