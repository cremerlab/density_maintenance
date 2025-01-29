#%%
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
exp_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
lit_prot = pd.read_csv('../../data/collated/collated_literature_total_protein.csv')
lit_rna = pd.read_csv('../../data/collated/collated_literature_total_RNA.csv')
fits = pd.read_csv('../../data/mcmc/rna_protein_per_cell_ppc_summary.csv')

#%%
fig, ax = plt.subplots(1, 2, figsize=(4, 2))

# Plot the literature data
for g, d in lit_prot.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt) 

for g, d in lit_rna.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    ax[1].plot(d['growth_rate_hr'], d['fg_rna_per_cell'], **fmt)

# Plot our data
fmt = size.viz.style_point('This Study')
ax[0].plot(exp_data['growth_rate_hr'], exp_data['fg_protein_per_cell'], **fmt)
ax[1].plot(exp_data['growth_rate_hr'], exp_data['fg_rna_per_cell'], **fmt)

# Populate a legend 
ax[0].legend()
ax[1].legend()

# Plot our fits
for i, q in enumerate(['prot_per_cell_ppc', 'rna_per_cell_ppc']):
    d = fits[fits['quantity']==q]
    ax[i].fill_between(d['growth_rate_hr'], d['sig2_lower'], d['sig2_upper'],
                       color=cor['pale_black'], alpha=0.5)
    ax[i].fill_between(d['growth_rate_hr'], d['sig1_lower'], d['sig1_upper'],
                       color=cor['light_black'], alpha=0.5)
    ax[i].plot(d['growth_rate_hr'], d['median'], '-', lw=1, color=cor['primary_black'])