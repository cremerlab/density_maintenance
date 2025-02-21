#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the many literature data sets
lit_ms = pd.read_csv('../../data/collated/merged_mass_spectrometry.csv')
lit_rp = pd.read_csv('../../data/collated/collated_literature_rna_to_protein.csv')
lit_size = pd.read_csv('../../data/collated/collated_literature_size_data.csv')
lit_prot = pd.read_csv('../../data/collated/collated_literature_total_protein.csv')
lit_rna = pd.read_csv('../../data/collated/collated_literature_total_RNA.csv')

# Load our data sets
our_size = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
our_size = our_size[our_size['strain']=='wildtype']
our_rp = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')


#%% 
# Assign cog superclasses
rev_mapper = {'metabolism': ['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'],
          'signaling': ['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'],
          'info': ['J', 'A', 'K', 'L', 'B'],
          'other': ['R', 'S', 'X']}
mapper = {_v:k for k, v in rev_mapper.items() for _v in v}
lit_ms['cog_group'] = [mapper[v] for v in lit_ms['cog_letter'].values]
grouped = lit_ms.groupby(['strain', 'source', 'growth_rate_hr', 'carbon_source', 'replicate', 'cog_group'])['mass_frac'].sum().reset_index()

#%%
# Set figure canvas for mass spec
fig, ax = plt.subplots(1, 4, figsize=(6, 1.75), sharex=True)
axes = {'metabolism': ax[0], 'info': ax[1], 'signaling': ax[2], 'other': ax[3]}
labels = ['COG: Metabolism', 'COG: Info. Storage & Processing', 'COG: Processes & Signaling', 'COG: Other']
# Plot data 
for g, d in grouped.groupby(['source', 'cog_group']):
    fmt = size.viz.style_point(g[0])
    _ax = axes[g[1]]
    _ax.plot(d['growth_rate_hr'], d['mass_frac'], **fmt)

# Add context
for i, a in enumerate(ax):
    a.set_xlim([0, 2.5])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_ylabel(f'{labels[i]}\nproteome mass fraction', fontsize=6)

ax[0].set_ylim([0.1, 0.7])
ax[1].set_ylim([0.1, 0.7])
ax[2].set_ylim([0, 0.3])
ax[3].set_ylim([0, 0.1])
ax[0].legend()
plt.tight_layout()
plt.savefig('./plots/figS2_MS_consistency.pdf', bbox_inches='tight')

#%% Set figure canvas for cell size
fig, ax = plt.subplots(1, 4, figsize=(6, 1.75))    
axes = {'length_um': ax[0], 'width_um': ax[1], 'volume_um3':ax[2], 'surface_area_um2': ax[3]}

for g, d in lit_size.groupby('source'):
    fmt = size.viz.style_point(g)
    for q, a in axes.items():
        a.plot(d['growth_rate_hr'], d[q], **fmt)
fmt = size.viz.style_point('This Study')
our_size.rename(columns={'volume_fL':'volume_um3'}, inplace=True)
for q, a in axes.items():
    a.plot(our_size['growth_rate_hr'], our_size[q], **fmt)

# Add context
for a in ax:
    a.set_xlim([0, 2.5])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylabel('average cell length [µm]', fontsize=6)
ax[0].set_ylim([1, 5.5])
ax[1].set_ylabel('average cell width [µm]', fontsize=6)
ax[1].set_ylim([0.4, 1.4])
ax[2].set_ylabel('average cell volume [µm$^{3}$]', fontsize=6)
ax[2].set_ylim([0, 7])
ax[3].set_ylabel('average cell surface area [µm$^{2}$]', fontsize=6)
ax[3].set_ylim([0, 22])
ax[0].legend()
plt.tight_layout()
plt.savefig('./plots/figS2_size_consistency.pdf', bbox_inches='tight')

#%% 
# Set figure canvas for bulk
fig, ax = plt.subplots(1, 3, figsize=(5, 1.75), sharex=True)

# Protein per cell
for g, d in lit_prot.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'
ax[0].plot(our_rp['growth_rate_hr'], our_rp['fg_protein_per_cell'], **fmt)

# RNA per cell
for g, d in lit_rna.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['fg_rna_per_cell'], **fmt)
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'
ax[1].plot(our_rp['growth_rate_hr'], our_rp['fg_rna_per_cell'], **fmt)

for g, d in lit_rp.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['rna_to_protein'], **fmt)
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'
ax[2].plot(our_rp['growth_rate_hr'], our_rp['fg_rna_per_cell']/our_rp['fg_protein_per_cell'], **fmt)

# Set context
for a in ax:
    a.set_xlim([0, 2.5])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.legend()
ax[0].set_ylabel('total RNA [fg / cell]', fontsize=6)
ax[0].set_ylim([0, 1000])
ax[1].set_ylabel('total protein [fg / cell]', fontsize=6)
ax[1].set_ylim([0, 500])
ax[2].set_ylabel('RNA-to-protein', fontsize=6)
ax[2].set_ylim([0, 0.7])

plt.tight_layout()
plt.savefig('./plots/figS2_bulk_consistency.pdf', bbox_inches='tight')
