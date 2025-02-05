#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import scipy.stats
import size.analytical
import size.viz 
import size.analytical
cor, pal = size.viz.matplotlib_style()

# Load data to motivate fixed aspect ratio
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']=='wildtype']

# Define the approximate aspect ratio and periplasmic width
ALPHA = 4
PERI_WIDTH = 0.025 # in um

# Define the range of cell widths over which to plot the size scaling
width_range =  np.linspace(0.4, 1.2, 100)
idx = np.where(width_range <= 0.8)[0][-1]
 
# Set up a figure showing aspect ratio and growth-rate dependent compartment scaling.
fig, ax = plt.subplots(1, 2, figsize=(4, 2))
ax = ax.ravel()

# Plot the aspect ratio and compute an empirical fit
fmt = size.viz.style_point('This Study')
ax[0].plot(data['growth_rate_hr'], data['length_um'] / data['width_um'], **fmt)
ax[0].hlines(ALPHA, 0, 2.5, linestyle='-', lw=1)

# Plot the cytoplasmic volume
total_vol = size.analytical.volume(ALPHA * width_range, width_range)
total_SA = size.analytical.surface_area(ALPHA * width_range, width_range)
peri_vol = size.analytical.envelope_volume(ALPHA * width_range, width_range, PERI_WIDTH)
cyto_vol = total_vol - peri_vol
mem_area = 2 * total_SA

ax[1].plot(width_range/width_range[idx], cyto_vol/cyto_vol[idx], '-', lw=1,
           color=cor['primary_black'])
ax[1].plot(width_range/width_range[idx], mem_area/mem_area[idx], '-', lw=1,
           color=cor['primary_blue'])
ax[1].plot(width_range/width_range[idx], peri_vol/peri_vol[idx], '--', lw=1,
           color=cor['primary_purple'])

# Add context
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\ell / w$\ncell aspecct ratio', fontsize=6)
ax[1].set_ylabel('$S(w)/S(w_0)$\nrelative compartment size', fontsize=6)
ax[1].set_xlabel('relative cell width\n$w / w_0$', fontsize=6)
ax[0].set_ylim([0, 5])
plt.tight_layout()
plt.savefig('./plots/figS2_size_scaling.pdf', bbox_inches='tight')
