# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Load the mass spec data
data = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')

cyto = data[data['envelope'] == False]
cyto = cyto.groupby(['dataset_name', 'condition',
                     'growth_rate_hr', 'ribosomal', 'metabolism'])['mass_frac'].sum().reset_index()
membrane = data[data['membrane'] == True]
membrane = membrane.groupby(['dataset_name', 'condition',
                             'growth_rate_hr'])['mass_frac'].sum().reset_index()
periplasm = data[data['periplasm'] == True]
periplasm = periplasm.groupby(['dataset_name', 'condition',
                               'growth_rate_hr'])['mass_frac'].sum().reset_index()

fig, ax = plt.subplots(1, 3, figsize=(6, 1.75))
for g, d in cyto[cyto['ribosomal'] == True].groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=mapper[g[0]]['c'],
               alpha=0.5)

for g, d in membrane.groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=mapper[g]['c'], alpha=0.5,
               label=g)

for g, d in periplasm.groupby(['dataset_name']):
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=mapper[g]['c'], alpha=0.5)
ax[1].legend()
ax[0].set_ylim([0, 0.5])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 0.13])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_title('allocation towards translation proteins', fontsize=6)
ax[1].set_title('allocation towards membrane proteins', fontsize=6)
ax[2].set_title('allocation towards periplasmic proteins', fontsize=6)
ax[0].set_ylabel('fraction of total proteome', fontsize=6)
ax[1].set_ylabel('fraction of total proteome', fontsize=6)
ax[2].set_ylabel('fraction of total proteome', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig1_allocation_trends.pdf')
