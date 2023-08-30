# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import glob

# Using information from matlab scripts, define the acceptable frame numbers
# to keep for valid analysis.
mu9_bounds = [-1.5, 3]
min_duration = 15
mats = [
    # ['Minimal Media 1ng.mat', {'time_bounds': [250, 450],
    #    'spike_time': 337.5833,
    #    'ignore_ids': [],
    #    'key': 1}],
    # ['Minimal Media 1ng.mat', {'time_bounds': [-100, 0],
    #    'ignore_ids': [],
    #    'key': 2,
    #    'spike_time': 337.5833}],
    # ['Minimal Media 2ng.mat', {'time_bounds': [400, 680],
    #    'ignore_ids': [541, 750, 831, 1493, 1884],
    #    'key': 1,
    #    'spike_time': 245.0667}],
    # ['Minimal Media 2ng.mat', {'time_bounds': [-150, 0],
    #    'ignore_ids': [541, 750, 831, 1493, 1884],
    #    'key': 2,
    #    'spike_time': 245.0667}],
    # ['Minimal Media ppGpp0.mat', {'time_bounds': [5, 1E6],
    #   'ignore_ids': [],
    #   'key': 1,
    #   'spike_time': 5}],
    ['Minimal Media 100uM.mat', {'time_bounds': [400, 1E6],
                                 'ignore_ids': [],
                                 'key': 1,
                                 'spike_time':512.033}],
    ['Minimal Media 100uM.mat', {'time_bounds': [-100, 0],
                                 'ignore_ids': [],
                                 'key': 2,
                                 'spike_time':512.033}]]


# 'Minimal Media ppGpp0.mat': [5, 1E6]}

# %%
cell_df = pd.DataFrame([])
agg_df = pd.DataFrame([])

for file in mats:
    bound = file[1]
    mat = scipy.io.loadmat(
        f'../../../data/literature/Buke2022/{file[0]}', squeeze_me=True, simplify_cells=True)
    conc = file[0].split(' ')[-1].split('.mat')[0]
    oe = 'none'
    inducer = 0
    if bound['key'] == 1:
        if 'ng' in conc:
            oe = 'relA'
            inducer = float(conc[:-2])
        if 'uM' in conc:
            oe = 'meshI'
            inducer = float(conc[:-2])
    for i, m in enumerate(mat['schnitzcells']):
        times = m['time'] - bound['spike_time']
        if (type(times) is float) | (type(times) is int):
            times = np.array(times)
        if ((times >= bound['time_bounds'][0]).all() & (times <= bound['time_bounds'][1]).all()) & (m['completeCycle'] != 0) & (m['interDivTime'] > min_duration) & (i+1 not in bound['ignore_ids']):
            avg_width = np.mean(m['rp_width'])
            _df = pd.DataFrame({'overexpression': oe,
                                'inducer_conc': inducer,
                                'key': bound['key'],
                                'width': avg_width,
                                'lam': np.mean(m['av_mu_rp']),
                                'file': file[0]},
                               index=[0])
            cell_df = pd.concat([cell_df, _df], sort=False)
agg_size = cell_df.groupby(
    ['overexpression', 'inducer_conc', 'key']).mean().reset_index()
agg_size['inducer_conc'] = agg_size['inducer_conc'].values.astype(float)
# %%
agg_size
# %%
# Load the ppGpp measurements and calculate the relative value to the uninduced sample average
ppGpp = pd.read_csv(
    '../../../data/literature/Buke2022/buke2022_ppGpp_concentrations.csv')
ppGpp_0 = ppGpp[(ppGpp['inducer_conc'] == 0) & (
    ppGpp['overexpression'] == 'meshI')]['ppGpp_per_biomass'].mean()
ppGpp = ppGpp.groupby(['inducer_conc', 'overexpression']).mean().reset_index()
ppGpp['rel_ppGpp'] = ppGpp['ppGpp_per_biomass'] / ppGpp_0

# Use the flux-parity theory to compute the ribosomal allocation.
phiRb_0 = 0.118  # Computed as 0.4558 * 0.26, from caption of Fig S3C
phiO = 0.55  # From flux parity theory, Chure 2023.
ppGpp['calc_phiRb'] = (
    (1 + ppGpp['rel_ppGpp'] * (1 - phiO - phiRb_0)/(phiRb_0))/(1 - phiO))**-1


# # Unify the dataframes
# for g, d in ppGpp.groupby(['overexpression', 'inducer_conc']):
#     agg_size.loc[(agg_size['overexpression'] == g[0]) & (agg_size['inducer_conc'] == g[1]),
#                  'calc_phiRb'] = d['calc_phiRb'].values[0]
#     agg_size.loc[(agg_size['overexpression'] == g[0]) & (agg_size['inducer_conc'] == g[1]),
#                  'calc_phiRb'] = d['calc_phiRb'].values[0]

agg_size['phiRb'] = [0.4558 * 0.26, 0.4558 * 0.31]
agg_size['growth_rate_hr'] = agg_size['lam'] * np.log(2)
agg_size['source'] = 'Buke et al. 2022'
agg_size.to_csv(
    '../../../data/literature/Buke2022/buke2022_processed.csv', index=False)

# %%
ppGpp
