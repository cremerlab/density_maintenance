#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import scipy.stats
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load literature data for RNA/protein and size
rp_data = pd.read_csv('../../data/collated/collated_literature_rna_to_protein.csv')
size_data = pd.read_csv('../../data/collated/collated_literature_size_data.csv') 

# Set the figure canvas
fig, ax = plt.subplots(1, 2, figsize=(4, 2))

# Plot the RP data
for g, d in rp_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['rna_to_protein'], **fmt)

# Plot the cell volume data 
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['volume_um3'], **fmt) 

# Compute empirical linear regressions 
lam_range = np.linspace(0, 2.5, 200)
rp_popt = scipy.stats.linregress(rp_data['growth_rate_hr'], 
                                 rp_data['rna_to_protein'])
vol_popt = scipy.stats.linregress(size_data['growth_rate_hr'],
                                  np.log(size_data['volume_um3']))
rp_fit = rp_popt[1] + rp_popt[0] * lam_range
vol_fit = np.exp(vol_popt[1] + vol_popt[0] * lam_range)

# Plot the fits
ax[0].plot(lam_range, rp_fit, 'k--', lw=1)
ax[1].plot(lam_range, vol_fit, 'k--', lw=1)

# Set context
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('RNA-to-protein', fontsize=6)
ax[1].set_ylabel('average cell volume [µm$^3$]', fontsize=6)
ax[0].set_ylim([0, 0.6])
ax[1].set_yscale('log')
ax[1].set_ylim([0.1, 10])
ax[0].legend()
ax[1].legend()
plt.savefig('./plots/figS1_growth_law_plots.pdf', bbox_inches='tight')