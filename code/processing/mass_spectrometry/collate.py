"""
This script collates all growth curves, cell size measurements, total protein,
and total RNA measurements into single dataframes.
"""
# %%
import scipy.stats
import numpy as np
import pandas as pd
import glob

# Load the total protein and total RNA data and merge
tot_prot = pd.read_csv('./2024-05-29_r2/raw/2024-05-29_r2_biuret.csv')
tot_RNA = pd.read_csv('./2024-05-30_r1/raw/2024-05-30_r1_total_RNA.csv')

# Apply transformations to total RNA to account for cell loss
tot_RNA['net_od_600nm'] = 1.5 * (tot_RNA['harvest_od_600nm'] -
                                 tot_RNA['residual_od_600nm']) - 1.2 * tot_RNA['wash_od600nm']
tot_RNA['od_260nm_per_od_600nm'] = tot_RNA['od_260nm'] / \
    tot_RNA['net_od_600nm']

# Transform total protein optical density measurements to account for
tot_prot['net_od_600nm'] = 1.5 * \
    (tot_prot['harvest_od_600nm'] - tot_prot['residual_od_600nm'])
tot_prot['od_555nm_per_od_600nm'] = tot_prot['od_555nm'] * \
    (0.4/1.5) / tot_prot['net_od_600nm']

# Merge the total protein and total RNA data
merged = pd.merge(tot_prot, tot_RNA, on=[
                  'strain', 'carbon_source', 'replicate', 'inducer_conc', 'date_collected'])
merged = merged[['strain', 'carbon_source', 'replicate', 'inducer_conc',
                 'date_collected', 'od_555nm_per_od_600nm', 'od_260nm_per_od_600nm']]
merged.to_csv('./collated_protein_RNA_measurements.csv', index=False)


# %%
# Create a dictionary mapping date to strain, carbon source, and replicate to the
mapper = {g[:-1]: g[-1] for (g, _) in merged.groupby(
    ['strain', 'carbon_source', 'inducer_conc', 'replicate', 'date_collected'])}

size_data = pd.concat([pd.read_csv(f)
                      for f in glob.glob('./*_r1/processed/*sizes.csv')])
growth_data = pd.concat([pd.read_csv(f) for f in glob.glob(
    './*_r1/processed/*growth_curves_processed.csv')])
growth_data.loc[growth_data['strain'].isnull(), 'strain'] = 'wildtype'
growth_data.loc[growth_data['inducer_conc'].isnull(), 'inducer_conc'] = 0
valid_sizes = pd.DataFrame()
valid_growth = pd.DataFrame()

for (k, v) in mapper.items():
    _valid_sizes = size_data[(size_data['strain'] == k[0]) & (size_data['carbon_source'] == k[1]) &
                             (size_data['replicate'] == k[3]) & (size_data['inducer_conc'] == k[2]) &
                             (size_data['date'] == v)]
    _valid_growth = growth_data[(growth_data['strain'] == k[0]) & (growth_data['carbon_source'] == k[1]) &
                                (growth_data['replicate'] == k[3]) & (growth_data['inducer_conc'] == k[2]) &
                                (growth_data['date'] == v)]
    valid_sizes = pd.concat([valid_sizes, _valid_sizes])
    valid_growth = pd.concat([valid_growth, _valid_growth])

valid_sizes = valid_sizes[['date', 'strain', 'carbon_source', 'replicate', 'inducer_conc',
                           'width_median', 'width_var', 'length', 'volume', 'surface_area', 'surface_to_volume',
                           'cell_id', 'image']]
valid_growth = valid_growth[['date', 'strain', 'carbon_source', 'replicate', 'inducer_conc',
                             'elapsed_time_hr', 'od_600nm']]
valid_sizes.to_csv('./collated_single_cell_size_measurements.csv', index=False)
valid_growth.to_csv('./collated_growth_curves.csv', index=False)

# %%
# Infer growth rates
growth_rates = pd.DataFrame([]) 
keys = ['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate']
for g, d in valid_growth.groupby(keys):
    d['elapsed_time_hr'] -= d['elapsed_time_hr'].min()
    popt = scipy.stats.linregress(d['elapsed_time_hr'].values,
                                  np.log(d['od_600nm']))
    _df_dict = {k:v for k, v in zip(keys, g)} 
    _df_dict['growth_rate_hr'] = popt[0]
    _df_dict['od_init'] = np.exp(popt[1])
    _df_dict['growth_rate_hr_std'] = popt[-2]
    growth_rates = pd.concat([growth_rates, pd.DataFrame(_df_dict, index=[0])])
growth_rates.to_csv('./collated_growth_rates.csv', index=False)

# %%
grouped_sizes = valid_sizes.groupby(['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate'])[
    ['width_median', 'width_var', 'length', 'volume', 'surface_area', 'surface_to_volume']].mean().reset_index()
grouped_sizes.to_csv('./aggregated_size_measurements.csv', index=False)
