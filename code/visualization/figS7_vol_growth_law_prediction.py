#%%
import numpy as np
import pandas as pd 
import scipy.stats
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load experimental dataset
exp_data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
wt_data = exp_data[exp_data['strain']=='wildtype']

# Restrict to wildtype for a simple exponential fit.
popt = scipy.stats.linregress(wt_data['growth_rate_hr'], np.log(wt_data['volume_fL']))

# Precict the volume based on the growth rate.
exp_data['predicted_volume_fL'] = np.exp(popt[1] + popt[0] * exp_data['growth_rate_hr'])
wt_data = exp_data[exp_data['strain']=='wildtype']

# Define colors for induction conditions
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'D', 'glucoseCAA': 's'}

# Set figure canvas
fig, ax = plt.subplots(1, 2, figsize=(4, 1.5))

# Define the wildtype points
fmt = size.viz.style_point('This Study')
fmt['markerfacecolor'] = cor['pale_black']
fmt['markeredgewidth'] = 0
fmt['alpha'] = 0.75
fmt['markersize'] = 5

# Plot the data and fit
ax[0].plot(wt_data['growth_rate_hr'], wt_data['volume_fL'], **fmt)

lam_range = np.linspace(0, 2.5, 50)
fit = np.exp(popt[1] + popt[0] * lam_range)
ax[0].plot(lam_range, fit, '-', color='black', lw=0.5)


# Plot wildtype data for predicted vs measured
ax[1].plot(wt_data['volume_fL'], wt_data['predicted_volume_fL'], **fmt)

# Plot induction conditions
for g, d in exp_data[exp_data['strain']!='wildtype'].groupby(['carbon_source', 'strain', 'inducer_conc']):
    marker = markers[g[0]] 
    markeredgecolor = cor[f'dark_{strain_colors[g[1]]}']
    color = cor[f'{inducer_concs[g[2]]}{strain_colors[g[1]]}']
    fmt = {'marker':marker,
           'markeredgecolor':markeredgecolor,
           'color':color, 
           'alpha':0.85,
           'markersize':4,
           'linestyle':'none',
           'alpha': 0.75}
    ax[0].plot(d['growth_rate_hr'], d['volume_fL'], **fmt)
    ax[1].plot(d['volume_fL'], d['predicted_volume_fL'], **fmt)

# Plot a 1:1 identity line
ax[1].plot([0, 5], [0, 5], 'k-', lw=0.5)

# Add plot labels and context
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('cell volume [µm$^3$]', fontsize=6)
ax[1].set_xlabel('measured cell volume [µm$^3$]', fontsize=6)
ax[1].set_ylabel('predicted cell volume [µm$^3$]', fontsize=6)

# Set scaling
ax[0].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_yscale('log')


# Adjust limits 
ticks = [0.5, 1, 2, 4]
lims = [0.4, 4.1]
for i, a in enumerate(ax):
    if i > 0:
        a.set_xlim(lims)
        a.set_xticks(ticks)
        a.set_xticklabels(ticks)
    a.set_ylim(lims)
    a.set_yticks(ticks)
    a.set_yticklabels(ticks)



plt.tight_layout()
plt.savefig('./plots/figS7_vol_growth_law_prediction.pdf', bbox_inches='tight')