# %%
import numpy as np
import pandas as pd
import scipy.stats
import tqdm
data = pd.read_csv(
    '../../../data/growth_curves/growth_measurements_processed.csv')
summarized = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['strain', 'carbon_source',  'overexpression', 'inducer', 'inducer_conc', 'temperature', 'run_idx'])):
    popt = scipy.stats.linregress(
        d['elapsed_time_hr'].values, np.log(d['od_600nm']))
    slope = np.round(popt[0], decimals=3)
    intercept = np.round(popt[1], decimals=3)
    _df = pd.DataFrame([np.array([g[0], g[1], g[2], g[3], g[4], g[5], g[6], slope, np.exp(intercept)])],
                       columns=['strain', 'carbon_source', 'overexpression', 'inducer',  'inducer_conc', 'temperature', 'run_idx', 'growth_rate_hr', 'od_init'])

    summarized = pd.concat([summarized, _df], sort=False)
summarized.loc[summarized['inducer'] == 'chloramphenicol', 'inducer'] = 'cm'
summarized['temperature'] = summarized['temperature'].astype(float)
summarized['inducer_conc'] = summarized['inducer_conc'].values.astype(float)

# %%
# Restrict to only wildtype, relA, and mesh for now.
summarized = summarized[summarized['strain'].isin(['wildtype', 'lpp14']) &
                        summarized['overexpression'].isin(['none', 'relA', 'meshI']) &
                        summarized['inducer_conc'].isin([0.0, 1.0, 2.0, 10.0, 100])]
# %%
summarized.to_csv(
    '../../../data/summaries/summarized_growth_measurements.csv', index=False)
