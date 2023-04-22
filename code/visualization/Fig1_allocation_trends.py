# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Load the mass spec data
data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
data_complete = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')
cyto_cogs = data_complete[data_complete['envelope'] == False].groupby(['dataset_name', 'growth_rate_hr', 'condition', 'cog_class'])[
    'mass_frac'].sum().reset_index()
cyto_cogs
env_cogs = data_complete[data_complete['envelope'] == True].groupby(['dataset_name', 'growth_rate_hr', 'condition', 'cog_class'])[
    'mass_frac'].sum().reset_index()
# %%
fig, ax = plt.subplots(2, 1, figsize=(3, 3.5), sharex=True)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

for a in ax:
    a.set_ylim([-0.05, 1])
ax[0].set_ylabel('composition of\ncytoplasm proteome', fontsize=6)
ax[1].set_ylabel('composition of\nenvelope proteome', fontsize=6)
cogcor = {'Not Assigned': cor['light_black'], 'cellular processes and signaling': cor['dark_red'],
          'information storage and processing': cor['primary_blue'],
          'metabolism': cor['primary_green'],
          'poorly characterized': cor['primary_black']}

for g, d in cyto_cogs.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    d['frac'] = d['mass_frac']/d['mass_frac'].sum()
    for _g, _d in d.groupby(['cog_class']):
        ax[0].plot(_d['growth_rate_hr'], _d['frac'], mapper[g[0]]['m'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4, color=cogcor[_g], alpha=0.5)

for g, d in env_cogs.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    d['frac'] = d['mass_frac']/d['mass_frac'].sum()
    for _g, _d in d.groupby(['cog_class']):
        ax[1].plot(_d['growth_rate_hr'], _d['frac'], mapper[g[0]]['m'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4, color=cogcor[_g], alpha=0.5)
plt.subplots_adjust(hspace=0.1)
plt.savefig('../../figures/Fig1_compartment_cog_classification.pdf',
            bbox_inches='tight')
# %%

inner_membrane = data[data['localization'] == 'inner membrane']
inner_membrane = inner_membrane.groupby(['dataset_name', 'condition',
                                         'growth_rate_hr'])['mass_frac'].sum().reset_index()
outer_membrane = data[data['localization'] == 'outer membrane']
outer_membrane = outer_membrane.groupby(['dataset_name', 'condition',
                                         'growth_rate_hr'])['mass_frac'].sum().reset_index()

periplasm = data[data['localization'] == 'periplasm']
periplasm = periplasm.groupby(['dataset_name', 'condition',
                               'growth_rate_hr'])['mass_frac'].sum().reset_index()

fig, ax = plt.subplots(1, 3, figsize=(6, 1.75))
for g, d in inner_membrane.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=cor['light_blue'],
               alpha=0.5)

for g, d in periplasm.groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=cor['primary_purple'], alpha=0.5,
               label=g)

for g, d in outer_membrane.groupby(['dataset_name',  'growth_rate_hr', 'condition']):
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=cor['primary_blue'], alpha=0.5)
# ax[1].legend()
for a in ax:
    a.set_yticks([0, 0.05, 0.10, 0.15])
    a.set_ylim([0, 0.15])

for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_title('allocation towards\ninner membrane proteins', fontsize=6)
ax[1].set_title('allocation towards\nperiplasmic proteins', fontsize=6)
ax[2].set_title('allocation towards\nouter membrane proteins', fontsize=6)
ax[0].set_ylabel('fraction of total proteome', fontsize=6)
ax[1].set_ylabel('fraction of total proteome', fontsize=6)
ax[2].set_ylabel('fraction of total proteome', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig1_allocation_trends.pdf')
