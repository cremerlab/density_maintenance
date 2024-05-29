#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()

DATE = "2024-05-28"
OD_BOUNDS = [0.04, 0.45]

data = pd.read_csv(f'./raw/{DATE}_r1_growth_curves_raw.csv')
data = data[(data['od_600nm'] >= OD_BOUNDS[0]) &
            (data['od_600nm'] <= OD_BOUNDS[1])]
#%%
# Calculate the time in minutes since the start of the experiment
df = pd.DataFrame([])
for g, d in data.groupby(['replicate', 'carbon_source', 'strain', 'inducer_conc']):
    d = d.copy()
    d['clock_time'] = pd.to_datetime(d['clock_time'], format='%H:%M')
    d['elapsed_time_hr'] = d['clock_time'] - d['clock_time'].min()
    d['elapsed_time_hr'] = d['elapsed_time_hr'].dt.total_seconds() / 3600
    df = pd.concat([df, d])

# Add important identifiers
df['temperature_C'] = 37
df.drop(columns=['clock_time'], inplace=True)
df.to_csv(f'./processed/{DATE}_r1_growth_curves_processed.csv', index=False)

#%%
# Compute the growth rates using standard scipy stats linregress
stats_df = pd.DataFrame([])
for g, d in df.groupby(['strain', 'inducer_conc', 'replicate', 'carbon_source']):
    popt = scipy.stats.linregress(d['elapsed_time_hr'].values, np.log(d['od_600nm'].values))  
    _df = pd.DataFrame({'date': DATE,
                        'strain': g[0],
                        'carbon_source': g[-1],
                        'replicate': g[2],
                        'inducer_conc': g[1],
                        'growth_rate_hr': popt[0],
                        'od_init': np.exp(popt[1]),
                        'growth_rate_std': popt[4]},
                        index=[0])
    stats_df = pd.concat([stats_df, _df]) 

stats_df.to_csv(f'./processed/{DATE}_r1_growth_rates.csv', index=False)

#%%
# Plot the growth rates
df['axis'] = df.groupby(['strain', 'carbon_source', 'inducer_conc', 'replicate']).ngroup()
n_cols = 3
n_rows = int(np.ceil((df['axis'].max()+1) / n_cols))
resid = (n_cols * n_rows) - (df['axis'].max() + 1)
fig, ax = plt.subplots(n_rows, n_cols, figsize=(6, 2 * n_rows))
ax = ax.ravel()

for i in range(1, resid+1):
    ax[-i].axis(False)

for a in ax:
    a.set_xlabel('elapsed time [hr]', fontsize=6)
    a.set_ylabel('optical density at 600nm [a.u.]', fontsize=6)
for g, d in df.groupby(['strain', 'carbon_source', 'inducer_conc', 'replicate', 'axis']):
    _fit = stats_df[(stats_df['strain'] == g[0]) &
                    (stats_df['carbon_source'] == g[1]) &
                    (stats_df['inducer_conc'] == g[2]) &
                    (stats_df['replicate']==g[3])]
    time_range = np.linspace(0, d['elapsed_time_hr'].max(), 100)
    fit = _fit['od_init'].values[0] * np.exp(_fit['growth_rate_hr'].values[0] * time_range)
    ax[g[4]].plot(d['elapsed_time_hr'], d['od_600nm'], 'o', ms=4, color=cor['primary_black'])
    ax[g[4]].plot(time_range, fit, '-', color=cor['primary_red'], lw=1)
    ax[g[4]].set_title(f"{g[0]}, {g[1]} Ind{g[2]} r{g[3]}: $\lambda$={_fit['growth_rate_hr'].values[0]:0.2f} inv.hr.", fontsize=6)
for a in ax:
    a.set_yscale('log')
plt.tight_layout()
plt.savefig(f'./viz/{DATE}_r1_growth_curves.png')