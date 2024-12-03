#%%
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

# Load the data
data = pd.read_csv('../../data/compiled_measurements.csv')
data = data[data['strain']=='wildtype']

# Load the parameter summarys
pars = pd.read_csv('../../data/mcmc/fig3_parameter_summaries.csv')

# Load the trend quantiles
quants = pd.read_csv('../../data/mcmc/fig3_fits.csv')

#%%
# Plot of the trends and the fits
fig, ax = plt.subplots(1, 2, figsize=(3, 1))

# Set postitional and stying information for the quantiles
axes = {'phi_mem_ppc': [ax[0], 'blue'],
        'phi_peri_ppc': [ax[1], 'purple']}
perc_colors = {'2sig':'pale_', '1sig':'light_'}

# Iterate through each quantity and plot
for g, d in quants[quants['quantity'].isin(axes)].groupby('quantity'):
    a = axes[g][0]
    a.plot(d['phi_rib'], d['median_val'], '-', lw=1, color=cor[axes[g][1]])
    for k, v in perc_colors.items():
        c = cor[f'{v}{axes[g][1]}'] 
        a.fill_between(d['phi_rib'], d[f'{k}_lower'], d[f'{k}_upper'],
                       color=c, alpha=0.5)

# Plot markers
fmt = size.viz.style_point('This Study')
fmt['markeredgecolor'] = cor['primary_blue']
ax[0].plot(data['phi_rib'], data['phi_mem'], **fmt)
fmt['markeredgecolor'] = cor['primary_purple']
ax[1].plot(data['phi_rib'], data['phi_peri'], **fmt)

# Add context
for a in ax:
    a.set_xlabel('ribosomal proteome\nfraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6) 
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)

# Adjust axes
ax[0].set_xlim([0.1, 0.3])
ax[1].set_xlim([0.1, 0.3])
ax[0].set_ylim([0, 0.20])
ax[1].set_ylim([0, 0.15])

