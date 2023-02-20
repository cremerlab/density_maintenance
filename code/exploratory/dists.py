# %%
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import size.viz
cor, pal = size.viz.matplotlib_style()
sizes = pd.read_csv(
    '../processing/microscopy/size_measurement/output/compiled_size_measurements.csv')
sizes = sizes[(sizes['strain'].isin(['wildtype'])) &
              (sizes['overexpression'] == 'none')]
# sizes = sizes[sizes['overexpression'].isin(['rbsB'])]
# sizes['aspect_ratio'] = sizes['length'] / sizes['width_median']
x = sizes.groupby(['strain', 'carbon_source', 'overexpression',
                  'inducer_conc', 'date', 'run_no']).mean().reset_index()
x['aspect_ratio'] = x['length'] / x['width_median']
x.groupby(['carbon_source', 'overexpression', 'inducer_conc'])[
    'length'].agg(('mean', 'sem'))
# %%
fig, ax = plt.subplots(2, 1, figsize=(8, 8))
for g, d in sizes.groupby(['strain', 'overexpression', 'inducer_conc']):
    if g[-1] == 0:
        ax[0].hist(d['width_median'], bins=50, alpha=0.5,
                   label=f'{g[0]}, {g[1]}, {g[-1]}', density=True)
    else:
        if g[1] != 'none':
            ax[1].hist(d['width_median'], bins=50, alpha=0.5,
                       label=f'{g[0]}, {g[1]}, {g[-1]}', density=True)


ax[0].legend()
ax[1].legend()


# %%

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for g, d in sizes.groupby(['strain', 'overexpression', 'inducer_conc']):
    if (g[-1] in [0, 100]) & (g[1] != 'none'):
        ax.hist(d['width_median'], bins=50, alpha=0.5,
                label=f'{g[0]}, {g[1]}, {g[-1]}', density=True)
ax.legend()

# %%
sns.boxplot(sizes, y='width_median', x='overexpression', hue='inducer_conc')


# %%
_data = sizes[(sizes['inducer_conc'].isin([0, 10, 30, 50, 100]))]
_data = _data[(_data['width_median'] >= 0.25) &
              (_data['width_median'] <= 1.25)]
_data = _data[((_data['overexpression'] == 'none') & (_data['inducer_conc'] == 0)) |
              ((_data['overexpression'] == 'malE') & (_data['inducer_conc'] > 0))]
sns.boxplot(_data, y='width_median', x='inducer_conc', hue='overexpression')
