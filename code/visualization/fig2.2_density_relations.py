#%%
import pandas as pd
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the protein per cell data
lit_prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
prot_data = pd.read_csv('../../data/bulk_protein_per_cell.csv')

# Load the size data
lit_size = pd.read_csv('../../data/literature/full_literature_size_data.csv')
data = pd.read_csv('../../data/compiled_measurements.csv')
data = data[data['strain']=='wildtype']
data['source'] = 'This Study'

# Load the MCMC empirical summaries
quants = pd.read_csv('../../data/mcmc/fig2_empirical_quantities_summary.csv')
fits = pd.read_csv('../../data/mcmc/fig2_ppc_summary.csv')


#%%
fig, ax = plt.subplots(1, 1, figsize=(1.75, 1)) 

# Plot the inference
fit = fits[fits['quantity']=='prot_per_cell_ppc']
ax.fill_between(fit['growth_rate_hr'], fit['lower'], fit['upper'], 
                color=cor['pale_black'], alpha=0.75)
ax.plot(fit['growth_rate_hr'], fit['mean_val'], 'k-', lw=1)

# Plot the literature data
for g, d in lit_prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'].values, d['fg_protein_per_cell'], **fmt)

# Plot our experimental data
fmt = size.viz.style_point('This Study')
fmt['color'] = cor['primary_black']
fmt['lw'] =0.75 
fmt['markerfacecolor'] = 'w'
fmt['capsize'] = 0
ax.errorbar(prot_data['mean_growth_rate_hr'], prot_data['fg_prot_per_cell'],
            xerr=prot_data['std_growth_rate_hr'], yerr=prot_data['err_fg_prot_per_cell'],
            **fmt)


# Add context
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('$M_{prot}$\nprotein per cell [fg/fL]', fontsize=6)
ax.set_ylim([0, 1000])
ax.set_xlim([0, 2.2])
# ax.legend()
plt.savefig('./plots/fig2.2_mprot_small.pdf', bbox_inches='tight')

#%% Plot the masses inferred from the mass spectrometry data
fig, ax = plt.subplots(3, 1, figsize=(2, 3), sharex=True)
mapper = {'M_cyto': [ax[0], cor['light_black'], cor['black']],
          'M_peri': [ax[1], cor['light_purple'], cor['purple']],
          'M_mem': [ax[2], cor['light_blue'], cor['blue']]}
for g, d in quants[quants['quantity'].isin(mapper.keys())].groupby(['quantity', 'source']):
    fmt = size.viz.style_point(g[1]) 
    fmt['markerfacecolor'] = mapper[g[0]][1]
    lw = 0.5
    if g[1] == 'This Study':
        fmt['markerfacecolor'] = 'w'
        fmt['markeredgecolor'] = mapper[g[0]][2]
        lw = 1
    mapper[g[0]][0].vlines(d['growth_rate_hr'], d['lower'], d['upper'], lw=lw,
              color=fmt['markeredgecolor'])
    mapper[g[0]][0].plot(d['growth_rate_hr'], d['mean_val'], **fmt)

# Add context
ax[0].set_ylim([10, 1000])
ax[1].set_ylim([-10, 50])
ax[2].set_ylim([-10, 125])
ax[0].set_ylabel('$M_{cyto}$ [fg / cell]', fontsize=6)
ax[1].set_ylabel('$M_{peri}$ [fg / cell]', fontsize=6)
ax[2].set_ylabel('$M_{mem}$ [fg / cell]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].legend()
plt.savefig('./plots/fig2.2_protein_loads.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1, 2, figsize=(2, 1))


# Plot the inference
vol_fit = fits[fits['quantity']=='volume_ppc']
sa_fit = fits[fits['quantity']=='surface_area_ppc']
ax[0].fill_between(vol_fit['growth_rate_hr'], vol_fit['lower'], vol_fit['upper'], 
                color=cor['pale_black'], alpha=0.75)
ax[0].plot(vol_fit['growth_rate_hr'], vol_fit['mean_val'], 'k-', lw=1)

ax[1].fill_between(sa_fit['growth_rate_hr'], sa_fit['lower'], sa_fit['upper'], 
                color=cor['pale_black'], alpha=0.75)
ax[1].plot(sa_fit['growth_rate_hr'], sa_fit['mean_val'], 'k-', lw=1)

# Plot the literature data
for g, d in lit_size.groupby('source'):
    fmt = size.viz.style_point(g) 
    ax[0].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)

# Plot our size measurements
fmt = size.viz.style_point('This Study')
ax[0].plot(data['growth_rate_hr'], data['volume_fL'], **fmt)
ax[1].plot(data['growth_rate_hr'], data['surface_area_um2'], **fmt)

# Add context
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_xlim(0, 2.5)
ax[0].set_ylabel('$V_{cell}$\ntotal volume [fL]', fontsize=6)
ax[1].set_ylabel('$S_A$\nsurface area [µm$^{2}$]', fontsize=6)
ax[0].legend(fontsize=6)
plt.savefig('./plots/fig2.2_compartment_size_fits.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(3, 1, figsize=(2, 3), sharex=True)
mapper = {'rho_cyto': [ax[0], cor['light_black'], cor['black']],
          'rho_peri': [ax[1], cor['light_purple'], cor['purple']],
          'sigma_mem': [ax[2], cor['light_blue'], cor['blue']]}
for g, d in quants[quants['quantity'].isin(mapper.keys())].groupby(['quantity', 'source']):
    fmt = size.viz.style_point(g[1]) 
    fmt['markerfacecolor'] = mapper[g[0]][1]
    lw = 0.5
    if g[1] == 'This Study':
        fmt['markerfacecolor'] = 'w'
        fmt['markeredgecolor'] = mapper[g[0]][2]
        lw = 1
    mapper[g[0]][0].vlines(d['growth_rate_hr'], d['lower'], d['upper'], lw=lw,
              color=fmt['markeredgecolor'])
    mapper[g[0]][0].plot(d['growth_rate_hr'], d['mean_val'], **fmt)

# Add context
ax[0].set_ylim([0, 600])
ax[1].set_ylim([0, 250])
ax[2].set_ylim([0, 8])
ax[0].set_ylabel(r'$\rho_{cyto}$ [fg / fL]', fontsize=6)
ax[1].set_ylabel(r'$\rho_{peri}$ [fg / fL]', fontsize=6)
ax[2].set_ylabel(r'$\sigma_{mem}$ [fg / µm$^{2}$]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.savefig('./plots/fig2.2_protein_densities.pdf', bbox_inches='tight')
