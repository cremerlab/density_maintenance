#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Define constants
BETA = 2.19 

# Load the data, parameter samples, and wild-type predictions
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')


#%%
# Define a graphical language for the perturbations
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'D', 'glucoseCAA': 's'}

# Plot the changes in the growth rate, ribosome content, and SAV
fig, axes = plt.subplots(5, 3, figsize=(3.5, 1.5))

# Set the columnsj
cols = {'psi_cyto':0, 'psi_mem':1,  'psi_peri':2}
rows = {('meshI', 100): 0, ('meshI', 0): 1, 
        ('relA', 0): 2, ('relA', 2): 3, ('relA', 4):4}

# Iterate through and plot the individual replicates
for g, d in data[data['strain']!='wildtype'].groupby(['strain', 'inducer_conc', 'carbon_source']):
    for q, c in cols.items():

        # Slice out the axis
        ax = axes[rows[g[:2]], c]

        # Style the point
        fmt = {'marker': markers[g[2]],
               'markersize': 4,
               'markeredgecolor': cor[f'dark_{strain_colors[g[0]]}'],
               'markerfacecolor': cor[f'{inducer_concs[g[1]]}{strain_colors[g[0]]}'],
               'linestyle': 'none',
               'alpha':0.75
               }
        ax.plot(d[q], [0, 1], **fmt)

# Adjust ranges and set context
for i in range(5):
    axes[i, 0].set_xlim(0.7, 0.9)
    axes[i, 1].set_xlim(0.05, 0.25)  
    axes[i, 2].set_xlim(0, 0.14)

axes[-1, 0].set_xlabel('cytoplasmic\nproteome partition\n$\psi_{cyto}$', fontsize=6)
axes[-1, 1].set_xlabel('membrane\nproteome partitiion\n$\psi_{mem}$', fontsize=6)
axes[-1, 2].set_xlabel('periplasm\nproteome partition\n$\psi_{peri}$', fontsize=6)

# Add specific ticking
axes[-1, 0].set_xticks([0.7,  0.8, 0.9])
axes[-1, 1].set_xticks([0.05, 0.15, 0.25])
axes[-1, 2].set_xticks([0, 0.075, 0.14])
for i in range(4):
    for j in range(3):
        axes[i, j].set_xticks([])   


for a in axes.ravel():
    a.set_yticks([])
    a.set_ylim([-3, 3])
    a.set_facecolor('w')
    a.spines['bottom'].set_color('black')
    a.spines['bottom'].set_linewidth(0.25)
    a.tick_params(axis='x', direction='out', width=0.25, length=1.5, 
                  colors='black')
plt.savefig('./plots/figSX_partition_induction_trends.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(2, 3, figsize=(5,3))

wt = data[data['strain']=='wildtype']
quants = ['psi_cyto', 'psi_mem', 'psi_peri']
for i, q in enumerate(quants):
    for j, var in enumerate(['phi_rib', 'growth_rate_hr']):
        ax[j, i].plot(wt[var], wt[q], 'o', markeredgewidth=0,
                  color=cor['pale_black'], alpha=0.5, ms=4) 

   
# Iterate through and plot the individual replicates
for g, d in data[data['strain']!='wildtype'].groupby(['strain', 'inducer_conc', 'carbon_source']):
    for q, c in cols.items():
        # Style the point
        fmt = {'marker': markers[g[2]],
               'markersize': 4,
               'markeredgecolor': cor[f'dark_{strain_colors[g[0]]}'],
               'markerfacecolor': cor[f'{inducer_concs[g[1]]}{strain_colors[g[0]]}'],
               'linestyle': 'none',
               'alpha':0.75
               }
        for j, var in enumerate(['phi_rib', 'growth_rate_hr']):
            ax[j, c].plot(d[var], d[q], **fmt)

# Add plot context 
for i in range(3):
    ax[0, i].set_xlim([0.1, 0.25])
    ax[1, i].set_xlim([0, 1.5])
    ax[0, i].set_xlabel('ribosome proteome allocation\n$\phi_{rib}$', fontsize=6)
    ax[1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
for i in range(2):
    ax[i, 0].set_ylim([0.7, 0.9])
    ax[i, 0].set_ylabel('$\psi_{cyto}$\ncytoplasm\nproteome allocation', fontsize=6)
    ax[i, 1].set_ylabel('$\psi_{mem}$\nmembrane\nproteome allocation', fontsize=6)
    ax[i, 2].set_ylabel('$\psi_{peri}$\nmembrane\nproteome allocation', fontsize=6)
    ax[i, 1].set_ylim([0, 0.25])
    ax[i, 2].set_ylim([0,  0.15])
plt.tight_layout()
plt.savefig('./plots/figSX_partition_induction_dependence.pdf', bbox_inches='tight')