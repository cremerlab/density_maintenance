# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
cor, pal = size.viz.matplotlib_style()

# %%
mass_spec_data = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')
mass_spec_data

# %%
# Envelope proteins
periplasmic = mass_spec_data[mass_spec_data['periplasm'] == True]
cog_periplasmic = periplasmic.groupby(['cog_letter', 'dataset_name',
                                       'growth_rate_hr', 'condition']).sum().reset_index()
relative_periplasm = pd.DataFrame([])
for g, d in cog_periplasmic.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    d['relative_mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    relative_periplasm = pd.concat([relative_periplasm, d], sort=False)
relative_periplasm.sort_values(
    by='growth_rate_hr', inplace=True, ascending=False)


cmap = sns.color_palette('tab20', n_colors=len(
    cog_periplasmic['cog_letter'].unique()))
cmap = {k: v for k, v in zip(cog_periplasmic['cog_letter'].unique(), cmap)}
# ignore dataset difference, look only at schmidt
# sel = relative_periplasm[relative_periplasm['dataset_name']
#  == 'Schmidt et al. 2016']
fig, ax = plt.subplots(1, 1, figsize=(4, 6))
labels = []
ylabels = []
for i, (g, d) in enumerate(relative_periplasm.groupby(['growth_rate_hr', 'condition', 'dataset_name'])):
    ylabels.append(f'{g[0]:0.3f}')
    d.sort_values(by='relative_mass_frac', inplace=True, ascending=False)
    tot_frac = 0
    for _g, _d in d.groupby('cog_letter', sort=False):
        if _g in labels:
            label = '__nolegend__'
        else:
            label = _g
            labels.append(_g)
        p = ax.barh(i, _d['relative_mass_frac'], 0.7,
                    label=label, left=tot_frac, color=cmap[_g])

        # ax.hlines(i, tot_frac, tot_frac + _d['relative_mass_frac'], lw=10, label=label,
        #   color=cmap[_g])
        tot_frac += _d['relative_mass_frac'].values[0]

ax.set_yticks(np.arange(len(ylabels)))
ax.set_yticklabels(ylabels)
ax.legend(bbox_to_anchor=(1.08, 1.0))
ax.set_ylabel('growth rate [hr$^{-1}$]')
ax.set_xlabel('fraction of total periplasmic protein')
plt.savefig('./periplasmic_space_cogs.pdf', bbox_inches='tight')


# %%
# Do a reclassification
periplasmic = mass_spec_data[mass_spec_data['periplasm'] == True]
periplasmic['category'] = 'everything else'
periplasmic.loc[periplasmic['cog_letter'].isin(
    ['G', 'E', 'F', 'H', 'P', 'I']), 'category'] = 'transport processes\n[G, E, F, H, P, I]'
periplasmic.loc[periplasmic['cog_letter'].isin(
    ['M', 'D']), 'category'] = 'division and envelope biogenesis\n[M, D]'
periplasmic.loc[periplasmic['cog_letter'].isin(
    ['N', 'T']), 'category'] = 'motility & signaling\n[N, T]'

grouped_periplasmic = periplasmic.groupby(['category', 'dataset_name',
                                           'growth_rate_hr', 'condition']).sum().reset_index()
relative_periplasm = pd.DataFrame([])
for g, d in grouped_periplasmic.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    d['relative_mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    relative_periplasm = pd.concat([relative_periplasm, d], sort=False)
relative_periplasm.sort_values(
    by='growth_rate_hr', inplace=True, ascending=False)
cmap = sns.color_palette('hls', n_colors=len(
    grouped_periplasmic['category'].unique()))
cmap = {k: v for k, v in zip(grouped_periplasmic['category'].unique(), cmap)}

fig, ax = plt.subplots(1, 1, figsize=(4, 6))
labels = []
ylabels = []
order = ['transport processes\n[G, E, F, H, P, I]', 'division and envelope biogenesis\n[M, D]',
         'motility & signaling\n[N, T]', 'everything else']
for i, (g, d) in enumerate(relative_periplasm.groupby(['growth_rate_hr', 'condition', 'dataset_name'])):
    ylabels.append(f'{g[0]:0.3f}')
    d.sort_values(by='relative_mass_frac', inplace=True, ascending=False)
    tot_frac = 0

    for o in order:
        if o in labels:
            label = '__nolegend__'
        else:
            label = o
            labels.append(o)
        _d = d[d['category'] == o]
        p = ax.barh(i, _d['relative_mass_frac'], 0.7,
                    label=label, left=tot_frac, color=cmap[o])

        # ax.hlines(i, tot_frac, tot_frac + _d['relative_mass_frac'], lw=10, label=label,
        #   color=cmap[_g])
        tot_frac += _d['relative_mass_frac'].values[0]

ax.set_yticks(np.arange(len(ylabels)))
ax.set_yticklabels(ylabels)
ax.legend(bbox_to_anchor=(1.08, 1.0))
ax.set_ylabel('growth rate [hr$^{-1}$]')
ax.set_xlabel('fraction of total periplasmic protein')
plt.savefig('./periplasmic_space_subcategorization.pdf', bbox_inches='tight')
# %%
periplasmic = mass_spec_data[mass_spec_data['periplasm'] == True]
periplasmic['category'] = 'everything else'
periplasmic.loc[periplasmic['go_terms'].str.contains(
    'GO:0051119'), 'category'] = 'sugar transport'
# periplasmic.loc[(periplasmic['go_terms'].str.contains(
#     'GO:0055052') |
#     periplasmic['go_terms'].str.contains('GO:0043190') |
#     (periplasmic['go_terms'].str.contains('GO:0043192'))), 'category'] = 'ABC transporter'

grouped_periplasmic = periplasmic.groupby(['category', 'dataset_name',
                                           'growth_rate_hr', 'condition']).sum().reset_index()
relative_periplasm = pd.DataFrame([])
for g, d in grouped_periplasmic.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    d['relative_mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    relative_periplasm = pd.concat([relative_periplasm, d], sort=False)
relative_periplasm.sort_values(
    by='growth_rate_hr', inplace=True, ascending=False)
cmap = sns.color_palette('hls', n_colors=len(
    grouped_periplasmic['category'].unique()))
cmap = {k: v for k, v in zip(grouped_periplasmic['category'].unique(), cmap)}

fig, ax = plt.subplots(1, 1, figsize=(4, 6))
labels = []
ylabels = []
order = relative_periplasm['category'].unique()
for i, (g, d) in enumerate(relative_periplasm.groupby(['growth_rate_hr', 'condition', 'dataset_name'])):
    ylabels.append(f'{g[0]:0.3f}')
    d.sort_values(by='relative_mass_frac', inplace=True, ascending=False)
    tot_frac = 0

    for o in order:
        if o in labels:
            label = '__nolegend__'
        else:
            label = o
            labels.append(o)
        _d = d[d['category'] == o]
        p = ax.barh(i, _d['relative_mass_frac'], 0.7,
                    label=label, left=tot_frac, color=cmap[o])

        # ax.hlines(i, tot_frac, tot_frac + _d['relative_mass_frac'], lw=10, label=label,
        #   color=cmap[_g])
        if len(_d) > 0:
            tot_frac += _d['relative_mass_frac'].values[0]

ax.set_yticks(np.arange(len(ylabels)))
ax.set_yticklabels(ylabels)
ax.legend(bbox_to_anchor=(1.08, 1.0))
ax.set_ylabel('growth rate [hr$^{-1}$]')
ax.set_xlabel('fraction of total periplasmic protein')
# plt.savefig('./periplasmic_space_ABC.pdf', bbox_inches='tight')

# %%

# %%

# %%

# %%

# %%

# %%

# %%
# Look at just abc transporters, whether they are periplamic or not
mass_spec_data['category'] = 'others'

mass_spec_data.loc[(mass_spec_data['go_terms'].str.contains(
    'GO:0055052') |
    mass_spec_data['go_terms'].str.contains('GO:0043190') |
    (mass_spec_data['go_terms'].str.contains('GO:0043192'))), 'category'] = 'ABC transporter'
abcs = mass_spec_data[mass_spec_data['category'] == 'ABC transporter']
