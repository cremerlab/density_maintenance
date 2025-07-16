#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()
W_PERI = 0.025

# Load our datasets
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
densities = pd.read_csv('../../data/mcmc/empirical_densities_summary.csv')
rp_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
rp_data = rp_data[rp_data['strain']=='wildtype']
data = data[data['strain']=='wildtype']

# Literature data
lit_ms_data = pd.read_csv('../../data/collated/compiled_literature_allocation_assignments_wide.csv')
lit_size_data = pd.read_csv('../../data/collated/collated_literature_size_data.csv')
lit_prot_data = pd.read_csv('../../data/collated/collated_literature_total_protein.csv')
lit_rna_data = pd.read_csv('../../data/collated/collated_literature_total_RNA.csv') 

# Fit trends to literature data
log_prot_popt = scipy.stats.linregress(lit_prot_data['growth_rate_hr'], np.log(lit_prot_data['fg_protein_per_cell']))
log_rna_popt = scipy.stats.linregress(lit_rna_data['growth_rate_hr'], np.log(lit_rna_data['fg_rna_per_cell']))
log_vol_popt = scipy.stats.linregress(lit_size_data['growth_rate_hr'], np.log(lit_size_data['volume_um3']))
log_sa_popt = scipy.stats.linregress(lit_size_data['growth_rate_hr'], np.log(lit_size_data['surface_area_um2']))

# Compute quantities for literature data
lit_ms_data['prot_tot'] = np.exp(log_prot_popt[1] + log_prot_popt[0]*lit_ms_data['growth_rate_hr'])
lit_ms_data['rna_tot'] = np.exp(log_rna_popt[1] + log_rna_popt[0]*lit_ms_data['growth_rate_hr'])
lit_ms_data['vol_um3'] = np.exp(log_vol_popt[1] + log_vol_popt[0]*lit_ms_data['growth_rate_hr'])
lit_ms_data['sa_um2'] = np.exp(log_sa_popt[1] + log_sa_popt[0]*lit_ms_data['growth_rate_hr'])
lit_ms_data['rho_cyto'] = (lit_ms_data['psi_cyto'] * lit_ms_data['prot_tot'] + lit_ms_data['rna_tot']) / (lit_ms_data['vol_um3'] - W_PERI * lit_ms_data['sa_um2'])
lit_ms_data['sigma_mem'] = lit_ms_data['psi_mem'] * lit_ms_data['prot_tot'] / (2 * lit_ms_data['sa_um2'])
lit_ms_data['rho_peri'] = (lit_ms_data['psi_peri'] * lit_ms_data['prot_tot']) / (lit_ms_data['sa_um2'] * W_PERI)
lit_ms_data['empirical_kappa'] = lit_ms_data['rho_cyto'] / lit_ms_data['sigma_mem']

#%% Plot the literature data and empirical fits. 
fig, ax = plt.subplots(1, 4, figsize=(7, 1.5))

# Plot the lit data
for g, d in lit_prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

for g, d in lit_rna_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['fg_rna_per_cell'], **fmt)
   
for g, d in lit_size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)
    ax[3].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)


# Plot our data   
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'
ax[0].plot(rp_data['growth_rate_hr'], rp_data['fg_protein_per_cell'], **fmt)
ax[1].plot(rp_data['growth_rate_hr'], rp_data['fg_rna_per_cell'], **fmt)
fmt['marker'] = 'o' 
ax[2].plot(data['growth_rate_hr'], data['volume_fL'], **fmt)
ax[3].plot(data['growth_rate_hr'], data['surface_area_um2'], **fmt)

# Overlay the fits
lam_range = np.linspace(0, 2.5, 100)
ax[0].plot(lam_range, np.exp(log_prot_popt[1] + log_prot_popt[0]*lam_range), 'k-', lw=1, label='literature fit')
ax[1].plot(lam_range, np.exp(log_rna_popt[1] + log_rna_popt[0]*lam_range), 'k-', lw=1, label='literature fit')              
ax[2].plot(lam_range, np.exp(log_vol_popt[1] + log_vol_popt[0]*lam_range), 'k-', lw=1, label='literature fit')
ax[3].plot(lam_range, np.exp(log_sa_popt[1] + log_sa_popt[0]*lam_range), 'k-', lw=1, label='literature fit')


# Add  context
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.legend()

# Add ylabels
ax[0].set_ylabel('protein per cell [fg]', fontsize=6)
ax[1].set_ylabel('RNA per cell [fg]', fontsize=6)
ax[2].set_ylabel('cell volume [µm$^{3}$]', fontsize=6)
ax[3].set_ylabel('cell surface area [µm$^{2}$]', fontsize=6)

# Set bounds and scaling
ax[0].set_ylim([10, 1200])
ax[1].set_ylim([1, 600])
ax[2].set_ylim([0.1, 8])
ax[3].set_ylim([0.8, 30])

plt.tight_layout()
plt.savefig('./plots/figSX_literature_densities_fits.pdf', bbox_inches='tight')

#%% Plot the densities
fig, ax = plt.subplots(2, 4, figsize=(8.5, 3.5)) 
for g, d in lit_ms_data.groupby('source'):
    fmt = size.viz.style_point(g)

    # Cytoplasmic partition and density
    ax[0, 0].plot(d['growth_rate_hr'], d['psi_cyto'], **fmt)
    ax[1, 0].plot(d['growth_rate_hr'], d['rho_cyto'], **fmt)

    # Membrane partition and density
    fmt['markeredgecolor'] = cor['blue']
    fmt['color'] = cor['pale_blue']
    ax[0, 1].plot(d['growth_rate_hr'], d['psi_mem'], **fmt)
    ax[1, 1].plot(d['growth_rate_hr'], d['sigma_mem'], **fmt)

    # Periplasmic partition and density
    fmt['markeredgecolor'] = cor['purple']
    fmt['color'] = cor['pale_purple']
    ax[1, 2].plot(d['growth_rate_hr'], d['rho_peri'], **fmt)
    ax[0, 2].plot(d['growth_rate_hr'], d['psi_peri'], **fmt)

    # Empirical kappa
    fmt['markeredgecolor'] = cor['primary_red']
    fmt['color'] = cor['pale_red']
    ax[1, 3].plot(d['growth_rate_hr'], d['empirical_kappa'], **fmt)

# Overlay our experimental measurements
for g, d in data[data['strain']=='wildtype'].groupby('strain'):
    fmt = size.viz.style_point('This Study')
    fmt['color'] = cor['pale_black']

    # Cytoplasmic partition
    ax[0, 0].plot(d['growth_rate_hr'], d['psi_cyto'], **fmt)

    # Membrane density
    fmt['markeredgecolor'] = cor['blue']
    fmt['color'] = cor['pale_blue']
    ax[0, 1].plot(d['growth_rate_hr'], d['psi_mem'], **fmt)

    # Periplasmic density
    fmt['markeredgecolor'] = cor['purple']
    fmt['color'] = cor['pale_purple'] 
    ax[0, 2].plot(d['growth_rate_hr'], d['psi_peri'], **fmt)

# Overlay our empirical densities
wt = densities[densities['strain']=='wildtype']
fmt = size.viz.style_point('This Study')
fmt['color'] = cor['pale_black']
for _, row in wt[wt['quantity']=='rho_cyt_tot'].iterrows():
    ax[1, 0].vlines(row['growth_rate_hr'], row['sig2_lower'], row['sig2_upper'], 
                    color=fmt['markeredgecolor'], lw=0.5)
    ax[1, 0].vlines(row['growth_rate_hr'], row['sig1_lower'], row['sig1_upper'], 
                    color=fmt['markeredgecolor'], lw=1)
    ax[1, 0].plot(row['growth_rate_hr'], row['mean'], **fmt)


for _, row in wt[wt['quantity']=='sigma_mem'].iterrows():
    fmt['markeredgecolor'] = cor['blue']
    fmt['color'] = cor['pale_blue']
    ax[1, 1].vlines(row['growth_rate_hr'], row['sig2_lower'], row['sig2_upper'], 
                    color=fmt['markeredgecolor'], lw=0.5)
    ax[1, 1].vlines(row['growth_rate_hr'], row['sig1_lower'], row['sig1_upper'], 
                    color=fmt['markeredgecolor'], lw=1)
    ax[1, 1].plot(row['growth_rate_hr'], row['mean'], **fmt)


for _, row in wt[wt['quantity']=='rho_peri'].iterrows():
    fmt['markeredgecolor'] = cor['purple']
    fmt['color'] = cor['pale_purple']
    ax[1, 2].vlines(row['growth_rate_hr'], row['sig2_lower'], row['sig2_upper'], 
                    color=fmt['markeredgecolor'], lw=0.5)
    ax[1, 2].vlines(row['growth_rate_hr'], row['sig1_lower'], row['sig1_upper'], 
                    color=fmt['markeredgecolor'], lw=1)
    ax[1, 2].plot(row['growth_rate_hr'], row['mean'], **fmt)

for _, row in wt[wt['quantity']=='empirical_kappa'].iterrows():
    fmt['markeredgecolor'] = cor['red']
    fmt['color'] = cor['pale_red']
    ax[1, 3].vlines(row['growth_rate_hr'], row['sig2_lower'], row['sig2_upper'], 
                    color=fmt['markeredgecolor'], lw=0.5)
    ax[1, 3].vlines(row['growth_rate_hr'], row['sig1_lower'], row['sig1_upper'], 
                    color=fmt['markeredgecolor'], lw=1)
    ax[1, 3].plot(row['growth_rate_hr'], row['mean'], **fmt)

# Adjust axis limits
ax[0, -1].axis(False)
ax[0, 0].set_ylim([0.5, 1])
ax[0, 1].set_ylim([0, 0.25])
ax[0, 2].set_ylim([0, 0.15])
ax[1, 0].set_ylim([0, 600])
ax[1, 1].set_ylim([0, 4])
ax[1, 2].set_ylim([0, 300])
ax[1, 3].set_ylim([0, 600])
# Add ylabels
ax[0, 0].set_ylabel('$\psi_{cyto}$\ncytoplasmic proteome partition', fontsize=6)
ax[0, 1].set_ylabel('$\psi_{mem}$\nmembrane proteome partition', fontsize=6)
ax[0, 2].set_ylabel('$\psi_{peri}$\nperiplasmic proteome partition', fontsize=6)
ax[1, 0].set_ylabel(r'$\rho_{cyto}$' + '\n[fg / µm$^3$]', fontsize=6)   
ax[1, 1].set_ylabel(r'$\sigma_{mem}$' + '\n[fg / µm$^2$]', fontsize=6)
ax[1, 2].set_ylabel(r'$\rho_{peri}$' + '\n[fg / µm$^3$]', fontsize=6)
ax[1, 3].set_ylabel('cytoplasmic-membrane\ndensity ratio [µm$^{-1}$]', fontsize=6)
for i in range(4):
    if i < 3:
        ax[0, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    ax[1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    
plt.subplots_adjust(hspace=0.3, wspace=0.45)
ax[0, 0].legend()
plt.savefig('./plots/figSX_literature_density_comparison.pdf', bbox_inches='tight')




