#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz
cor, pal = size.viz.matplotlib_style()

# Define constants
BETA = 1/0.4558

# Load the data, parameter samples, and wild-type predictions
data = pd.read_csv('../../data/mcmc/predicted_SAV_summary.csv')

# Define a graphical language for the perturbations
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'D', 'glucoseCAA': 's'}

# Plot the changes in the growth rate, ribosome content, and SAV
fig, axes = plt.subplots(5, 3, figsize=(3.5, 1.5))

# Set the columnsj
cols = {'phi_rib':0, 'growth_rate_hr':1,  'measured_SAV':2}
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
    axes[i, 0].set_xlim(0.1, 0.2)
    axes[i, 1].set_xlim(0.3, 1.5)   
    axes[i, 2].set_xlim(5, 7.5)

axes[-1, 0].set_xlabel('ribosomal proteome allocation\n$\phi_{rib}$', fontsize=6)
axes[-1, 1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
axes[-1, 2].set_xlabel('surface-to-volume [µm$^{-1}$]\n$S_A/V$', fontsize=6)

# Add specific ticking
axes[-1, 0].set_xticks([0.1, 0.15, 0.20])
axes[-1, 1].set_xticks([0.3, 0.9,  1.5])
axes[-1, 2].set_xticks([5, 6, 7.5])
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
plt.savefig('./plots/fig4_induction_trends.pdf', bbox_inches='tight')

#%%  Plot the predicted vs measured curves
fig, ax = plt.subplots(1, 2, figsize=(4, 1.5))
for a in ax:
    a.plot([0, 12], [0, 12], 'k-', lw=0.5)

# Plot the wildtype data
wt_data = data[data['strain']=='wildtype']
for a in ax:
    fmt = size.viz.style_point('This Study')
    fmt['markerfacecolor'] = cor['pale_black']
    fmt['markeredgecolor'] = cor['pale_black']
    fmt['alpha'] = 0.5
    fmt['markersize'] = 4
    a.vlines(wt_data['measured_SAV'], wt_data['predicted_SAV_sig2_lower'], 
             wt_data['predicted_SAV_sig2_upper'], lw=0.5,
              color=cor['pale_black'], alpha=0.5)
    a.vlines(wt_data['measured_SAV'], wt_data['predicted_SAV_sig1_lower'], 
             wt_data['predicted_SAV_sig1_upper'], lw=1,
              color=cor['pale_black'], alpha=0.5)
    a.plot(wt_data['measured_SAV'], wt_data['predicted_SAV_mean'], **fmt)

rel = data[data['strain']=='relA']
for g, d in rel.groupby(['carbon_source', 'inducer_conc']):
    fmt = {'marker':markers[g[0]],
           'markeredgecolor': cor[f'dark_red'],
           'markeredgewidth': 0.5,
           'markersize': 3, 
           'markerfacecolor': cor[f'{inducer_concs[g[1]]}red'],
           'linestyle':'none', 
           'alpha': 0.75}
    ax[1].vlines(d['measured_SAV'], d['predicted_SAV_sig2_lower'], 
                 d['predicted_SAV_sig2_upper'], lw=0.5,
              color=cor['dark_red'])
    ax[1].vlines(d['measured_SAV'], d['predicted_SAV_sig1_lower'], 
                 d['predicted_SAV_sig1_upper'], lw=1,
              color=cor['dark_red'])

    ax[1].plot(d['measured_SAV'], d['predicted_SAV_mean'], **fmt)

rel = data[data['strain']=='meshI']
for g, d in rel.groupby(['carbon_source', 'inducer_conc']):
    fmt = {'marker':markers[g[0]],
           'markeredgecolor': cor[f'dark_gold'],
           'markeredgewidth': 0.5,
           'markersize': 3, 
           'markerfacecolor': cor[f'{inducer_concs[g[1]]}gold'],
           'linestyle':'none',
           'alpha': 0.75}
    ax[0].vlines(d['measured_SAV'], d['predicted_SAV_sig2_lower'], 
                 d['predicted_SAV_sig2_upper'], lw=0.5,
              color=cor['dark_gold'])
    ax[0].vlines(d['measured_SAV'], d['predicted_SAV_sig1_lower'], 
                 d['predicted_SAV_sig1_upper'], lw=1,
              color=cor['dark_gold'])

    ax[0].plot(d['measured_SAV'], d['predicted_SAV_mean'], **fmt)

# Add a shaded panel for the non-physiological regime
min_width = 0.55 # in µm
aspect_ratio = 4 # Assume, backed up with data.
max_SAV = ((12 * aspect_ratio) / (3 * aspect_ratio - 1)) / min_width

ax[1].fill_between([0, 8], max_SAV, 20, color=cor['pale_black'], alpha=0.5)
# Add context
for a in ax:
    a.set_xlim([3, 8])
    a.set_xlabel('measured surface-to-volume [µm$^{-1}$]\n$(S_A/V)_{meas}$', fontsize=6)
    a.set_ylabel('$(S_A/V)_{pred}$\npredicted\nsurface-to-volume [µm$^{-1}$]', fontsize=6)

ax[0].set_ylim([3, 8])
ax[1].set_ylim([3, 13])
plt.subplots_adjust(wspace=0.4)
plt.savefig('./plots/fig4_pred_vs_meas.pdf', bbox_inches='tight')