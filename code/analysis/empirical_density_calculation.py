# %%
import numpy as np
import pandas as pd
import scipy.stats
ms_data = pd.read_csv('../../data/literature/summarized_mass_fractions.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
# prot_data = prot_data[prot_data['source'] != 'Richa']
# prot_data = prot_data[prot_data['source'] != 'Wright & Lockhart 1964']
# prot_data = prot_data[(prot_data['source'] == 'Dennis & Bremer 1974') | (
# prot_data['source'] == 'Chohji et al. 1976')]
# prot_data = prot_data[]

# %%
sav_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['surface_to_volume'])
vol_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], np.log(size_data['volume_um3']))
prot_popt = scipy.stats.linregress(
    prot_data['growth_rate_hr'], np.log(prot_data['fg_protein_per_cell']))
width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um']
)
length_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['length_um']
)

# Compute the fits
lam_range = np.linspace(0, 3, 100)
sav_fit = sav_popt[1] + sav_popt[0] * lam_range
prot_fit = np.exp(prot_popt[1] + prot_popt[0] * lam_range)
vol_fit = np.exp(vol_popt[1] + vol_popt[0] * lam_range)
length_fit = length_popt[1] + length_popt[0] * lam_range
width_fit = width_popt[1] + width_popt[0] * lam_range
df = pd.DataFrame(np.array([lam_range, sav_fit, vol_fit, prot_fit, width_fit, length_fit]).T,
                  columns=['growth_rate_hr', 'surface_to_volume',
                           'volume', 'fg_protein_per_cell', 'width', 'length'])
df.to_csv('../../data/empirical_literature_trends.csv', index=False)

prot_data['volume'] = np.exp(
    vol_popt[1] + vol_popt[0] * prot_data['growth_rate_hr'])
prot_data['density'] = prot_data['fg_protein_per_cell'] / prot_data['volume']
prot_data.to_csv('../../data/literature/collated_total_protein_density.csv')
# %%
# Compute the densities
ms_data['volume'] = np.exp(vol_popt[1] + vol_popt[0]
                           * ms_data['growth_rate_hr'])
ms_data['surface_to_volume'] = sav_popt[1] + \
    sav_popt[0] * ms_data['growth_rate_hr']
ms_data['surface_area'] = ms_data['surface_to_volume'] * ms_data['volume']
ms_data['total_protein'] = np.exp(prot_popt[1] +
                                  prot_popt[0] * ms_data['growth_rate_hr'])
ms_data['width'] = width_popt[1] + width_popt[0] * ms_data['growth_rate_hr']
ms_data['mass_fg'] = ms_data['mass_frac'] * ms_data['total_protein']
ms_data.to_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv', index=False)
