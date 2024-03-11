# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./temperature_effects.stan')
phi_data = pd.read_csv(
    '../../../data/summaries/merged_rna_protein_measurements.csv')
phi_data = phi_data[
                    (phi_data['overexpression'] == 'none') &
                    (phi_data['inducer_conc'] == 0) &
                    (phi_data['temperature'] != 37)]  

growth = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')
growth = growth[(growth['overexpression'] == 'none') & 
                (growth['inducer_conc'] == 0) &
                (growth['temperature'] != 37) & 
                (growth['carbon_source'] != 'glucose')] 

size_data = pd.read_csv(
    '../../../data/summaries/summarized_size_measurements.csv')
size_data = size_data[(size_data['overexpression'] == 'none') &
                      (size_data['inducer_conc'] == 0) & 
                      (size_data['temperature'] != 37) & 
                      (size_data['carbon_source'] != 'glucose')]
mapper = {'glycerol': {25:1, 30:2, 42:3},
          'glucoseCAA': {25:4, 30:5, 42:6}}
_mapper = [['glycerol', 25], ['glycerol', 30], ['glycerol', 42],
           ['glucoseCAA', 25], ['gucoseCAA', 30], ['glucoseCAA', 42]]
for k, v in mapper.items():
    for _k, _v in v.items():
        for d in [phi_data, growth, size_data]:
            d.loc[(d['carbon_source']==k) & (d['temperature']==_k), 'idx'] = _v
# %%
size_vals = size_data.groupby('idx')[
    ['width_median', 'length', 'volume', 'surface_area', 'aspect_ratio']].agg((np.mean, np.std))
growth_vals = growth.groupby('idx')['growth_rate_hr'].agg((np.mean, np.std))
phi_vals = phi_data.groupby('idx')[
    ['ug_prot_per_biomass', 'ug_rna_per_biomass', 'phiRb']].agg((np.mean, np.std))

# %%
data_dict = {'J_cond': 6,
             'N_growth': len(growth),
             'growth_idx': growth['idx'].values.astype(int),
             'growth_rates': growth['growth_rate_hr'].values,
             'N_phi': len(phi_data),
             'phi_idx': phi_data['idx'].values.astype(int),

             'prot_per_biomass': phi_data['ug_prot_per_biomass'].values,
             'rna_per_biomass': phi_data['ug_rna_per_biomass'].values,
             'phiRb': phi_data['phiRb'].values,

             'N_size': len(size_data),
             'size_idx': size_data['idx'].values.astype(int),
             'widths': size_data['width_median'].values,
             'lengths': size_data['length'].values,
             'volumes': size_data['volume'].values,
             'surface_areas': size_data['surface_area'].values,
             'surface_to_volume': size_data['surface_to_volume'].values,
             'aspect_ratios': size_data['aspect_ratio'].values,
             }

_samples = model.sample(data=data_dict, adapt_delta=0.999)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(6, 1, figsize=(3, 6), sharex=True)
# ax[0].set_xlim([0, 8])
# Plot PPC distributions for non-size parameters
units = ['[µm]', '[µm]', '', '[µm$^{2}$]', '[µm$^3$]', 'µm$^{-1}$']
labels = [f'{k}\n{_k}' for k, v in mapper.items() for _k, _ in v.items()]
ticks = np.arange(0, 6) + 1

for i, q in enumerate(['width', 'length', 'aspect_ratio', 'surface_area', 'volume', 'surface_to_volume']):
    ppc = samples.posterior[f'{q}_ppc'].to_dataframe().reset_index()
    for g, d in ppc.groupby(f'{q}_ppc_dim_0'):
        ax[i].plot(1 + np.ones_like(d[::3]) * g, d[f'{q}_ppc'].values[::3],
                   '_', markeredgecolor=cor['primary_black'], markeredgewidth=0.1,
                   alpha=0.75)
    ax[i].set_title(q, fontsize=6)
    ax[i].set_ylabel(f"{q.replace('_', ' ')} {units[i]}", fontsize=6)
for g, d in size_data.groupby(['idx']):
    for i, v in enumerate(['width_median', 'length', 'aspect_ratio', 'surface_area', 'volume', 'surface_to_volume']):
        ax[i].plot(np.ones_like(d) * g, d[v], '_',
                   markeredgecolor=cor['primary_red'], markeredgewidth=0.5)
ax[4].set_xticks(ticks)
ax[4].set_xticklabels(labels, fontsize=6)

# %%
fig, ax = plt.subplots(4, 1, figsize=(2, 5), sharex=True)
# ax[0].set_xlim([0, 3])
# Plot PPC distributions for non-size parameters
for i, q in enumerate(['prot', 'rna', 'phiRb', 'growth_rates']):
    if q != 'growth_rates':
        if q == 'phiRb':
            n = q
        else:
            n = f'{q}_per_biomass'
    else:
        n = 'growth_rates'
    if q == 'cells':
        mod = 1E9
    else:
        mod = 1
    ppc = samples.posterior[f'{n}_ppc'].to_dataframe().reset_index()

    for g, d in ppc.groupby(f'{n}_ppc_dim_0'):
        ax[i].plot(1 + np.ones_like(d[::3]) * g, d[f'{n}_ppc'].values[::3] / mod,
                   '_', markeredgecolor=cor['primary_black'], markeredgewidth=0.1,
                   alpha=0.75)

ax[0].plot(phi_data['idx'],
           phi_data['ug_prot_per_biomass'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[0].set_title('total protein', fontsize=6)
ax[0].set_ylabel('prot\n[µg / OD$_{600nm}$]', fontsize=6)


ax[1].plot(phi_data['idx'],
           phi_data['ug_rna_per_biomass'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[1].set_title('total RNA', fontsize=6)
ax[1].set_ylabel('RNA\n[µg / OD$_{600nm}$]', fontsize=6)

ax[2].plot(phi_data['idx'],
           phi_data['phiRb'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[2].set_title('$\phi_{Rb}$', fontsize=6)
ax[2].set_ylabel('ribosomal allocation', fontsize=6)


ax[3].plot(growth['idx'],
           growth['growth_rate_hr'], '_',
           markeredgecolor=cor['primary_red'], markeredgewidth=1)
ax[3].set_title('growth rates', fontsize=6)
ax[3].set_ylabel('growth rate\n[hr$^{-1}$]', fontsize=6)
ax[3].set_xticks(ticks)
ax[3].set_xticklabels(labels)
# ax[2].set_xticks([1, 2])
# ax[2].set_xticklabels(['relA, 1', 'meshI, 100'], fontsize=6)

# %%
pars = ['width', 'length', 'volume', 'surface_area',
        'surface_to_volume', 'aspect_ratio', 'growth_rate', 'phiRb']
units = ['µm', 'µm', 'fL', 'µm^2', 'µm^-1', 'none', 'hr^-1', 'none']
par_df = pd.DataFrame([])
post_df = pd.DataFrame([])
for i, p in enumerate(pars):
    post = samples.posterior[f'{p}_mu'].to_dataframe().reset_index()
    for g, d in post.groupby(f'{p}_mu_dim_0'):
        med_val = np.median(d[f'{p}_mu'].values)
        mean_val = np.mean(d[f'{p}_mu'].values)
        percs = np.percentile(d[f'{p}_mu'].values, (2.5, 97.5))
        c, t = _mapper[g]
        _df = pd.DataFrame({
            'strain': 'wildtype',
            'carbon_source': c,
            'temperature': t,   
            'overexpression': 'none',
            'inducer': 0,
            'quantity': p,
            'median_value': med_val,
            'mean_value': mean_val,
            '2.5%': percs[0],
            '97.5%': percs[1],
            'unit': units[i]},
            index=[0])
        par_df = pd.concat([par_df, _df], sort=False)
        _post_df = d.copy()
        _post_df['strain'] = 'wildtype'
        _post_df['overexpression'] = 'none'
        _post_df['inducer_conc'] = 0
        _post_df['carbon_source'] = c
        _post_df['temperature'] = t
        _post_df['quantity'] = p
        _post_df['unit'] = units[i]
        _post_df.rename(columns={f'{p}_mu': 'value'}, inplace=True)
        _post_df.drop(columns=[f'{p}_mu_dim_0'], inplace=True)
        post_df = pd.concat([_post_df, post_df], sort=False)

par_df.to_csv(
    '../../../data/mcmc/temperature_perturbation_parameter_summaries.csv', index=False)
post_df.to_csv(
    '../../../data/mcmc/temperature_perturbation_samples.csv', index=False)
# %%
