#%%
import numpy as np 
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the protein per cell data
lit_prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
prot_data = pd.read_csv('../../data/bulk_protein_per_cell.csv')

# Perform an exponential fit for demonstration
lam = np.concatenate([lit_prot_data['growth_rate_hr'].values, prot_data['mean_growth_rate_hr'].values]).flatten()
m_prot = np.concatenate([lit_prot_data['fg_protein_per_cell'].values, prot_data['fg_prot_per_cell'].values])
popt = scipy.stats.linregress(lam, np.log(m_prot))

# Compute the fit range
lam_range = np.linspace(0, 2.2, 100)
fit = np.exp(popt[1] + popt[0] * lam_range)

#%%
fig, ax = plt.subplots(1, 1, figsize=(1.75, 1)) 
for g, d in lit_prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'].values, d['fg_protein_per_cell'], **fmt)

fmt = size.viz.style_point('This Study')
fmt['color'] = cor['primary_black']
fmt['lw'] =0.75 
fmt['markerfacecolor'] = 'w'
fmt['capsize'] = 0
ax.errorbar(prot_data['mean_growth_rate_hr'], prot_data['fg_prot_per_cell'],
            xerr=prot_data['std_growth_rate_hr'], yerr=prot_data['err_fg_prot_per_cell'],
            **fmt)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('$M_{prot}$\nprotein per cell [fg/fL]', fontsize=6)
ax.plot(lam_range, fit, 'k--', lw=1)
ax.set_ylim([0, 1000])
ax.set_xlim([0, 2.2])
# ax.legend()
plt.savefig('./plots/fig2.2_mprot_small.pdf', bbox_inches='tight')

#%% Plot the masses inferred from the mass spectrometry data



#%%
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']
