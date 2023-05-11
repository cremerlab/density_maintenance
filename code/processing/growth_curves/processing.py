# %%
import pandas as pd

# Load the dataset
data = pd.read_csv(
    '../../../data/growth_curves/growth_measurements_raw.csv')

# Filter on OD bounds
data = data[data['od_600nm'] <= 0.6]

# convert clock time to elapsed time
data['clock_time'] = pd.to_datetime(data['clock_time'].values)
processed = []
for g, d in data.groupby(['run_idx']):
    _d = d.copy()
    _d['elapsed_time_hr'] = (
        d['clock_time'] - d['clock_time'].min()).dt.total_seconds() / 3600
    if len(_d) >= 3:
        processed.append(_d)
processed = pd.concat(processed, sort=False)
processed.to_csv(
    '../../../data/growth_curves/growth_measurements_processed.csv', index=False)
