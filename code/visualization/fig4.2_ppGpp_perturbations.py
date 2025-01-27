#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz
cor, pal = size.viz.matplotlib_style()
# Define constants
BETA = 1/0.4558

# Load the data, parameter samples, and wild-type predictions
data = pd.read_csv('../../data/mcmc/perturbation_predicted_SAV_summary.csv')
data



# par_samples = pd.read_csv('../../data/mcmc/fig3_inference_samples.csv')
# quants = pd.read_csv('../../data/mcmc/fig3_fits.csv')
#%%
# Extract dependent quantities
phi_rib = data['phi_rib']
phi_mem = data['phi_mem']
phi_peri = data['phi_peri']

# Storage vectors for means and bounds
mean_vals = []
sig2_upper = [] 
sig2_lower = [] 
sig1_upper = [] 
sig1_lower = []

# Define the percentiles to compute
percs = [(2.5, 97.5), (16, 84)]

# Iterate through each measurement
for i in range(len(data)):
    # Compute the distribution of predicted SAV
    sav_pred = phi_mem[i] * par_samples['kappa'] / (2 * (1 + BETA * phi_rib[i] - phi_mem[i] - phi_peri[i]))

    # Compute the percentiles
    mean_vals.append(sav_pred.mean())
    for p, (low, high) in zip(percs, [(sig2_lower, sig2_upper), 
                                      (sig1_lower, sig1_upper)]):
        _p = np.percentile(sav_pred, p)     
        low.append(_p[0])
        high.append(_p[1])

# Update the dataframe
data['pred_sav_mean'] = mean_vals
data['1sig_lower'] = sig1_lower
data['1sig_upper'] = sig1_upper
data['2sig_lower'] = sig2_lower
data['2sig_upper'] = sig2_upper

#%% Define a graphical language for the perturbations
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'D', 'glucoseCAA': 's'}

#%% Plot the changes in the growth rate, ribosome content, and SAV
fig, axes = plt.subplots(5, 3, figsize=(3.5, 1.5))

# Set the columnsj
cols = {'growth_rate_hr':0, 'phi_rib':1, 'sav_inv_um':2}
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
    axes[i, 0].set_xlim(0.3, 1.5)
    axes[i, 1].set_xlim(0.1, 0.2)
    axes[i, 2].set_xlim(5, 7.5)

axes[-1, 0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
axes[-1, 1].set_xlabel('ribosomal\nproteome fraction $\phi_{rib}$', fontsize=6)
axes[-1, 2].set_xlabel('surface-to-volume [µm$^{-1}$]\n$S_A/V$', fontsize=6)

# Add specific ticking
axes[-1, 0].set_xticks([0.3, 0.9,  1.5])
axes[-1, 1].set_xticks([0.1, 0.15, 0.20])
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
fig, ax = plt.subplots(1, 2, figsize=(4, 1.25))
for a in ax:
    a.plot([0, 12], [0, 12], 'k--', lw=0.5)

# Plot the wildtype data
wt_data = data[data['strain']=='wildtype']
for a in ax:
    fmt = size.viz.style_point('This Study')
    fmt['markerfacecolor'] = cor['pale_black']
    fmt['markeredgecolor'] = cor['pale_black']
    fmt['alpha'] = 0.5
    fmt['markersize'] = 3
    a.vlines(wt_data['sav_inv_um'], wt_data['2sig_lower'], wt_data['2sig_upper'], lw=0.5,
              color=cor['pale_black'], alpha=0.5)
    a.plot(wt_data['sav_inv_um'], wt_data['pred_sav_mean'], **fmt)

rel = data[data['strain']=='relA']
for g, d in rel.groupby(['carbon_source', 'inducer_conc']):
    fmt = {'marker':markers[g[0]],
           'markeredgecolor': cor[f'dark_red'],
           'markeredgewidth': 0.5,
           'markersize': 3, 
           'markerfacecolor': cor[f'{inducer_concs[g[1]]}red'],
           'linestyle':'none'}
    ax[0].vlines(d['sav_inv_um'], d['2sig_lower'], d['2sig_upper'], lw=0.5,
              color=cor['dark_red'],)
    ax[0].plot(d['sav_inv_um'], d['pred_sav_mean'], **fmt)

mesh = data[data['strain']=='meshI']
for g, d in mesh.groupby(['carbon_source', 'inducer_conc']):
    fmt = {'marker':markers[g[0]],
           'markeredgecolor': cor[f'dark_gold'],
           'markeredgewidth': 0.5,
           'markersize': 3,
           'markerfacecolor': cor[f'{inducer_concs[g[1]]}gold'],
           'linestyle':'none'}
    ax[1].vlines(d['sav_inv_um'], d['2sig_lower'], d['2sig_upper'], lw=0.5,
              color=cor['dark_gold'],)
    ax[1].plot(d['sav_inv_um'], d['pred_sav_mean'], **fmt)

# Add a shaded panel for the non-physiological regime
min_width = 0.55 # in µm
aspect_ratio = 4 # Assume, backed up with data.
max_SAV = ((12 * aspect_ratio) / (3 * aspect_ratio - 1)) / min_width

for a in ax:
    a.fill_betweenx([0, 20], max_SAV, 10, color=cor['pale_black'])
# Add context
for a in ax:
    a.set_xlim([3, 9])
    a.set_xlabel('measured surface-to-volume [µm$^{-1}$]\n$(S_A/V)_{meas}$', fontsize=6)
    a.set_ylabel('$(S_A/V)_{pred}$\npredicted\nsurface-to-volume [µm$^{-1}$]', fontsize=6)

ax[0].set_ylim([3, 13])
ax[1].set_ylim([3, 8])
plt.savefig('./plots/fig4_pred_vs_meas.pdf', bbox_inches='tight')