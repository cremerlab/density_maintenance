# %%
import numpy as np
import pandas as pd
import scipy.stats
import tqdm
data = pd.read_csv(
    '../../../data/growth_curves/growth_measurements_processed.csv')

summarized = pd.DataFrame([])

for g, d in tqdm.tqdm(data.groupby(['strain', 'carbon_source', 'run_no', 'date'])):
    popt = scipy.stats.linregress(
        d['elapsed_time_hr'].values, np.log(d['od_600nm']))
    slope = popt[0]
    intercept = popt[1]
    _df = pd.DataFrame([np.array([g[0], g[1], g[2], g[3], slope, np.exp(intercept)])],
                       columns=['strain', 'carbon_source', 'run_no', 'date', 'growth_rate_hr', 'od_init'])
    summarized = pd.concat([summarized, _df], sort=False)

# %%
summarized.to_csv(
    '../../../data/summaries/summarized_growth_measurements.csv', index=False)