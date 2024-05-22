#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()

DATE = "2024-05-17"
OD_BOUNDS = [0.04, 0.45]

data = pd.read_csv('./raw/2024-05-17_r1_growth_curves_raw.csv')
data = data[(data[''od_600nm''] >= OD_BOUNDS[0]) &
            (data[''od_600nm''] <= OD_BOUNDS[1])]
# Calculate the time in minutes since the start of the experiment
df = pd.DataFrame([])
for g, d in data.groupby(['replicate', 'carbon_source']):
    d = d.copy()
    d['clock_time'] = pd.to_datetime(d['clock_time'], format='%H:%M')
    d['elapsed_time_hr'] = d['clock_time'] - d['clock_time'].min()
    d['elapsed_time_hr'] = d['elapsed_time_hr'].dt.total_seconds() / 3600
    # d['elapsed_time_hr'] = d['elapsed_time_hr'].values.astype(float)
    df = pd.concat([df, d])

# Add important identifiers
df['overexpression'] = None
df['inducer'] = None
df['inducer_conc'] = 0
df['temperature_C'] = 37
df.drop(columns=['clock_time'], inplace=True)
df.to_csv('./processed/2024-05-17_r1_growth_curves_processed.csv', index=False)

#%%
# Compute the growth rates using standard scipy stats linregress
stats_df = pd.DataFrame([])
OD_BOUNDS = [0.04, 0.45]
df = df[(df['od_600nm'] >= OD_BOUNDS[0]) &
        (df['od_600nm'] <= OD_BOUNDS[1])]

for g, d in df.groupby(['date', 'carbon_source', 'replicate']):
    popt = scipy.stats.linregress(d['elapsed_time_hr'].values, np.log(d['od_600nm'].values))  
    _df = pd.DataFrame({'date': g[0],
                        'carbon_source': g[1],
                        'replicate': g[2],
                        'growth_rate_hr': popt[0],
                        'od_init': np.exp(popt[1]),
                        'growth_rate_std': popt[4]},
                        index=[0])
    stats_df = pd.concat([stats_df, _df]) 

stats_df.to_csv('./processed/2024-05-17_r1_growth_rates.csv', index=False)

#%%
# Plot the growth rates
fig, ax = plt.subplots(1, 2, figsize=(4, 2), sharex=True)
for a in ax:
    a.set_xlabel('elapsed time [hr]', fontsize=6)
    a.set_ylabel('optical density at 600nm [a.u.]', fontsize=6)

for g, d in df.groupby(['carbon_source', 'replicate']):
    _fit = stats_df[(stats_df['carbon_source'] == g[0]) &
                    (stats_df['replicate'] == g[1])]
    time_range = np.linspace(0, d['elapsed_time_hr'].max(), 100)
    fit = _fit['od_init'].values[0] * np.exp(_fit['growth_rate_hr'].values[0] * time_range)
    ax[g[1]-1].plot(d['elapsed_time_hr'], d['od_600nm'], 'o', ms=4, color=cor['primary_black'])
    ax[g[1]-1].plot(time_range, fit, '-', color=cor['primary_red'], lw=1)
    ax[g[1]-1].set_title(f"{g[0]}, r{g[1]}: $\lambda$={_fit['growth_rate_hr'].values[0]:0.2f} inv.hr.", fontsize=6)
for a in ax:
    a.set_yscale('log')
plt.tight_layout()
plt.savefig('./viz/2024-05-17_r1_growth_curves.png')