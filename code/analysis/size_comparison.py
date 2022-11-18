# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
sizes = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')


glucose = sizes[sizes['carbon_source'] == 'LB']


for g, d in glucose.groupby(['date']):
    plt.hist(d['width_median'], bins=100, alpha=0.5, label=g)
plt.legend()

# %%

_sizes = sizes[sizes['date'] ==
               '2022-09-29_'].groupby(['carbon_source']).mean().reset_index()
__sizes = sizes[sizes['date'] !=
                '2022-09-29_'].groupby(['carbon_source']).mean().reset_index()

fig, ax = plt.subplots(1, 1)
labs = []
ticks = []
i = 0
for g, d in __sizes.groupby(['carbon_source']):
    try:
        _d = _sizes[_sizes['carbon_source'] == g]
        ax.plot(i, _d['width_median'], 'X')
        print(
            f"difference for {g} is {0 - _d['width_median'].values[0]} µm")
    except:
        a = 1
    ax.plot(i, d['width_median'], 'o')
    labs.append(g)
    ticks.append(i)
    i += 1

ax.set_xticks(ticks)
ax.set_xticklabels(labs)
