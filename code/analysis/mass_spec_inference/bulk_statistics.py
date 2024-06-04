# %%
import numpy as np
import pandas as pd
import scipy.stats

# Load the biuret data and calculate the calibration curve data
biuret_cal = pd.read_csv(
    '../../processing/mass_spectrometry/2024-05-28_r2/raw/2024-05-28_biuret_calibration.csv')
biuret_slope, biuret_inter = scipy.stats.linregress(biuret_cal['od_555nm'],
                                                    biuret_cal['nominal_BSA_conc_ug_mL'])[:2]

# Load the total protein and RNA data
rp_data = pd.read_csv(
    '../../processing/mass_spectrometry/collated_protein_RNA_measurements.csv')
rp_data.rename(columns={'date_collected': 'date'}, inplace=True)

# Calculate the RNA per od
rp_data['ug_RNA_per_biomass'] = rp_data['od_260nm_per_od_600nm'] * 31
rp_data['ug_prot_per_biomass'] = rp_data['od_555nm_per_od_600nm'] * \
    biuret_slope + biuret_inter
rp_data['RNA_to_protein'] = rp_data['ug_RNA_per_biomass'] / \
    rp_data['ug_prot_per_biomass']
rp_data['phi_Rb'] = 0.4558 * rp_data['RNA_to_protein']
rp_data = rp_data[['date', 'strain', 'carbon_source',
                   'inducer_conc', 'replicate', 'RNA_to_protein', 'phi_Rb']]
# %%
# Load the growth data, fit the exponential growth curve, and report the growth rate
growth_data = pd.read_csv(
    '../../processing/mass_spectrometry/collated_growth_curves.csv')
growth_rates = pd.DataFrame([])

for g, d in growth_data.groupby(['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate']):
    popt = scipy.stats.linregress(d['elapsed_time_hr'],  np.log(d['od_600nm']))
    _df = pd.DataFrame({'date': g[0],
                        'strain': g[1],
                        'carbon_source': g[2],
                        'inducer_conc': g[3],
                        'replicate': g[4],
                        'growth_rate_hr': popt[0],
                        'growth_rate_hr_std': popt[4]},
                       index=[0])
    growth_rates = pd.concat([growth_rates, _df], sort=False)

# %% Merge the growth rate, RNA/protein, and size data
size_data = pd.read_csv(
    '../../processing/mass_spectrometry/aggregated_size_measurements.csv')
merged = pd.merge(rp_data, growth_rates, on=[
                  'date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged = pd.merge(merged, size_data, on=[
                  'date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])
merged.rename({'width_median': 'mean_width',
               'length': 'mean_length',
               'volume': 'mean_volume',
               'surface_area': 'mean_surface_area',
               'surface_to_volume': 'mean_surface_to_volume'}, axis=1, inplace=True)
# %%
wt = merged[merged['strain'] == 'wildtype']
plt.plot(wt['growth_rate_hr'], wt['RNA_to_protein'], 'o')
