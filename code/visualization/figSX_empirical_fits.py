#%%
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
exp_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
lit_prot = pd.read_csv('../../data/collated/collated_literature_total_protein.csv')
lit_rna = pd.read_csv('../../data/collated/collated_literature_total_RNA.csv')

#%%
fig, ax = plt.subplots(1, 2, figsize=(4, 2))

# Plot the literature data
for g, d in lit_prot.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt) 
for g, d in lit_rna.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['fg_rna_per_cell'], **fmt)

# Plot our data
fmt = size.viz.style_point('This Study')
ax[0].plot(exp_data['growth_rate_hr'], exp_data['fg_protein_per_cell'], **fmt)
ax[1].plot(exp_data['growth_rate_hr'], exp_data['fg_rna_per_cell'], **fmt)
    
