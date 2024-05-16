#%%
import numpy as np 
import pandas as pd 
import scipy.stats
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('./raw/2024-05-17_spectrophotometer_calibration.csv')

# Generate a summation figure
fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.set_xscale('log')
ax.set_yscale('log')
# Add labels
ax.set_xlabel('dilution factor', fontsize=6)
ax.set_ylabel('$OD_{600nm}$\noptical density [a.u.]', fontsize=6)
for g, d in data.groupby('vessel_type'):
    ax.plot(d['dilution_factor'], d['od600nm'], 'o', label=g, 
            alpha=0.75)
ax.legend()

plt.savefig('./output/2024-05-17_od_v_dilution.png')

#%%
# Restrict the data to only the linear regime
# data = data[(data['od600nm'] >= 0.04) & (data['od600nm'] <= 0.4)]
# For each dilution factor, compute the OD relative to the cuvette
df = pd.DataFrame([])
for g, d in data.groupby('dilution_factor'):
    if '1cm cuvette' not in d['vessel_type'].values:
        break
    cuvette = d[d['vessel_type'] == '1cm cuvette']['od600nm'].values[0]     
    for k in ['16mm tube', '20mm tube', '25mm tube']:
        _d = d[d['vessel_type'] == k]
        if len(_d) == 0:
            break
        tube = _d['od600nm'].values[0]
        _df = pd.DataFrame({'dilution_factor': g, 
                        'vessel_type': k, 
                        'tube_od600nm': tube,
                        'cuvette_od600nm': cuvette,
                        'relative_od600nm': tube/cuvette,
                        }, 
                        index=[0])
        df = pd.concat([df, _df])


#%%
# Perform a linear regression
fit_df = pd.DataFrame([])
for g, d in df.groupby('vessel_type'):
    d = d[(d['cuvette_od600nm'] >= 0.04) & (d['cuvette_od600nm'] <= 0.45) &
          (d['tube_od600nm'] >= 0.04) & (d['tube_od600nm'] <= 0.45)]
    popt = scipy.stats.linregress(d['tube_od600nm'], d['cuvette_od600nm'])
    slope = popt[0]
    err = popt[-1]
    _df = pd.DataFrame({'vessel_type': g,
                        'cal_factor': slope,
                        'cal_factor_err': err,
                        'cal_intercept': popt[1]},
                        index = [0])
    fit_df = pd.concat([fit_df, _df])
# Save the result
fit_df.to_csv('./output/2024-05-17_calibration_factors.csv', index=False)

#%% 
# Create a figure showing the fit
fig, ax = plt.subplots(1, 3, figsize=(6, 2), sharex=True)
axes = {'16mm tube': ax[0], '20mm tube': ax[1], '25mm tube': ax[2]}
tube_range = np.linspace(0, 2, 100)
for g, d in df.groupby('vessel_type'):

    # Compute the fit
    factor = fit_df[fit_df['vessel_type'] == g]['cal_factor'].values[0]
    intercept = fit_df[fit_df['vessel_type'] == g]['cal_intercept'].values[0]
    fit = intercept + factor*tube_range

    # Subselect only the points used in the fitting    
    sel_points = d[(d['cuvette_od600nm'] >= 0.04) & (d['cuvette_od600nm'] <= 0.45) &
                   (d['tube_od600nm'] >= 0.04) & (d['tube_od600nm'] <= 0.45)]
    axes[g].plot(d['tube_od600nm'], d['cuvette_od600nm'], 'o', label='test points')
    axes[g].plot(sel_points['tube_od600nm'], sel_points['cuvette_od600nm'], 'o', 
                 color='white', markeredgecolor=cor['primary_red'],
                 markeredgewidth=1.25, ms=5, zorder=1000, label='fit points')
    axes[g].plot(tube_range, fit, '-', color=cor['primary_red'], lw=1.5, label='fit')
    
    # Add the title and factor
    axes[g].set_title(f'{g}: {factor:.2f}', fontsize=6)

for a in ax:
    a.set_xlabel('tube OD$_{600nm}$', fontsize=6)  
    a.set_ylabel('cuvette OD$_{600nm}$', fontsize=6)
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_ylim([1E-2, 2])
ax[0].legend()
plt.tight_layout()
plt.savefig('./output/2024-05-17_calibration_fit.png')
