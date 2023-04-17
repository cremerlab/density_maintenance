# %%
import numpy as np
import pandas as pd
import scipy.stats

ms_data = pd.read_csv('../../data/literature/summarized_mass_fractions.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')

# %%
sav_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['surface_to_volume'])
vol_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], np.log(size_data['volume_um3']))
prot_popt = scipy.stats.linregress(
    prot_data['growth_rate_hr'], prot_data['fg_protein_per_cell'])

# Compute the fits
lam_range = np.linspace(0, 3, 100)
sav_fit = sav_popt[1] + sav_popt[0] * lam_range
prot_fit = prot_popt[1] + prot_popt[0] * lam_range
vol_fit = np.exp(vol_popt[1] + vol_popt[0] * lam_range)
df = pd.DataFrame(np.array([lam_range, sav_fit, vol_fit, prot_fit]).T, 
                  columns=['growth_rate_hr', 'surface_to_volume', 
                           'volume', 'fg_protein_per_cell'])
df.to_csv('../../data/empirical_literature_trends.csv', sort=False)

#%%
# Compute the densities
ms_data['volume'] = np.exp(vol_popt[1] + vol_popt[0] * ms_data['growth_rate_hr'])
ms_data['surface_to_volume'] = sav_popt[1] + sav_popt[0] * ms_data['growth_rate_hr']
ms_data['total_protein'] = prot_popt[1] + prot_popt[0] * ms_data['growth_rate_hr']
ms_data['mass_fg'] = ms_data['mass_frac'] * ms_data['total_protein']
ms_data.to_csv('../../data/literature/collated_mass_fractions_empirics.csv', index=False)