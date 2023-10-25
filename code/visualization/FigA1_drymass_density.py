#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import scipy.stats
import size.viz 
cor, pal = size.viz.matplotlib_style()
ms_data = pd.read_csv('../../data/literature/collated_mass_fractions_empirics.csv')
drymass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')
prot = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
dna = pd.read_csv('../../data/literature/collated_dna_protein_ratio.csv')

# Do an empirical exponential fit to the protein data
popt = scipy.stats.linregress(prot['growth_rate_hr'], np.log(prot['fg_protein_per_cell']))

ms_data.head()

#%%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for g, d in dna.groupby('source'):
    fmt = size. 
    ax.plot()