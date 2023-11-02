# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
import tqdm
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./SAV_model_inference.stan')

size_data = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
size_data = size_data[(size_data['source'] != 'Si et al. 2017') &
                      (size_data['source'] != 'Taheri-Araghi et al. 2015')]
prot_data = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
mem_data = ms_data[ms_data['localization'] == 'membrane']
peri_data = ms_data[ms_data['localization'] == 'periplasm']

phiRb_data = pd.read_csv(
    '../../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']

density_data = pd.read_csv(
    '../../../data/literature/collated_drymass_densities.csv')


N_pred = 300
phiRb_range = np.linspace(0.05, 0.45, N_pred)
lam_range = np.linspace(0, 2.5, N_pred)
data_dict = {
    'N_pred': N_pred,
    'pred_phiRb_range': phiRb_range,
    'pred_lam_range': lam_range,

    'N_size': len(size_data),
    'size_lam': size_data['growth_rate_hr'].values,
    'surface_areas': size_data['surface_area_um2'].values,
    'aspect_ratios': size_data['length_um'].values / size_data['width_um'].values,

    'N_prot': len(prot_data),
    'prot_lam': prot_data['growth_rate_hr'].values,
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,

    'N_drymass': len(density_data),
    'drymass_lam': density_data['growth_rate_hr'].values,
    'drymass_density': density_data['drymass_density_fg_fL'].values,

    'N_ms': len(mem_data),
    'ms_lam': mem_data['growth_rate_hr'].values,
    'phi_mem': mem_data['mass_frac'].values,
    'phi_peri':  peri_data['mass_frac'].values,

    'N_phiRb': len(phiRb_data),
    'phiRb': phiRb_data['mass_fraction'].values,
    'phiRb_lam': phiRb_data['growth_rate_hr'].values,
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)


# %%
# Summaries the parameter percentiles
pars = ['rho_mem_mu', 'phi_mem_mu', 'm_peri',
        'drymass_mu', 'alpha_mu', 'kappa', 'log_prot_intercept', 'log_prot_slope']
full_post = samples.posterior[pars].to_dataframe().reset_index()
full_post['idx'] = 1
full_post[pars].to_csv(
    '../../../data/mcmc/theory_parameter_posterior_samples.csv', index=False)

# Compute the parameter percentiles
lowers = np.array([2.5, 12.5, 62.5, 50])
uppers = 100 - lowers
labels = ['95%', '75%', '25%', 'median']
perc = size.viz.compute_percentiles(
    full_post, pars, 'idx', lower_bounds=lowers, upper_bounds=uppers, interval_labels=labels)

perc[['upper', 'lower', 'quantity', 'interval']].to_csv(
    '../../../data/mcmc/theory_parameter_summaries_longform.csv', index=False)

wide_df = pd.DataFrame([])
for g, d in perc.groupby(['quantity']):
    med = d[d['interval'] == 'median']['lower'].values[0]
    perc = d[d['interval'] == '95%']
    _df = pd.DataFrame({'quantity': g,
                        'median_value': med,
                        '97.5%': perc['lower'].values[0],
                        '2.5%': perc['upper'].values[0]}, index=[0])
    wide_df = pd.concat([wide_df, _df], sort=False)
wide_df.to_csv(
    '../../../data/mcmc/theory_parameter_summaries_wide.csv', index=False)

# %%
# Summarize the predicted growth-rate dependencies
pars = ['pred_lam_prot', 'volume_pred', 'length_pred', 'aspect_ratio_pred', 'width_pred',
        'm_peri_pred', 'rho_peri_pred', 'rho_mem_pred', 'rho_dry_pred', 'phi_mem_pred', 'phi_Rb_pred', 'phi_peri_lam_pred',
        'SA_fit']
lam_par_percs = pd.DataFrame([])
for i, p in enumerate(tqdm.tqdm(pars)):
    post = samples.posterior[p].to_dataframe().reset_index()
    for j, lam in enumerate(lam_range):
        post.loc[post[f'{p}_dim_0'] == j, 'growth_rate_hr'] = lam
    percs = size.viz.compute_percentiles(post, p, 'growth_rate_hr',
                                         lower_bounds=lowers,
                                         upper_bounds=uppers,
                                         interval_labels=labels)
    lam_par_percs = pd.concat([lam_par_percs, percs], sort=False)
lam_par_percs.to_csv(
    '../../../data/mcmc/theory_growth_rate_prediction_summaries.csv', index=False)

# %%
# Summarize theory prediction
pars = ['SAV_theory', 'SAV_theory_simple', 'width_theory', 'width_theory_simple', 'phi_peri_pred', 'phi_mem_pred']
theo_par_percs = pd.DataFrame([])
for i, p in enumerate(tqdm.tqdm(pars)):
    post = samples.posterior[p].to_dataframe().reset_index()
    for j, phiRb in enumerate(phiRb_range):
        post.loc[post[f'{p}_dim_0'] == j, 'phiRb'] = phiRb
    percs = size.viz.compute_percentiles(post, p, 'phiRb',
                                         lower_bounds=lowers,
                                         upper_bounds=uppers,
                                         interval_labels=labels)
    theo_par_percs = pd.concat([theo_par_percs, percs], sort=False)
theo_par_percs.to_csv(
    '../../../data/mcmc/theory_phiRb_prediction_summaries.csv', index=False)

# %%
# Summarize the empirical caluclations for the mass spectrometry data
pars = ['ms_m_peri', 'ms_rho_peri', 'ms_rho_mem', 'ms_kappa']
ms_par_percs = pd.DataFrame([])
for i, p in enumerate(tqdm.tqdm(pars)):
    post = samples.posterior[p].to_dataframe().reset_index()
    for j in range(len(mem_data)):
        post.loc[post[f'{p}_dim_0'] == j,
                 'source'] = mem_data['dataset_name'].values[j]
        post.loc[post[f'{p}_dim_0'] == j,
                 'growth_rate_hr'] = mem_data['growth_rate_hr'].values[j]
    percs = size.viz.compute_percentiles(post, p, ['source', 'growth_rate_hr'],
                                         lower_bounds=lowers,
                                         upper_bounds=uppers,
                                         interval_labels=labels)
    ms_par_percs = pd.concat([ms_par_percs, percs], sort=False)
ms_par_percs.to_csv(
    '../../../data/mcmc/mass_spec_empirical_summaries_longform.csv', index=False)
# %%
ms_par_wide = pd.DataFrame([])
for g, d in ms_par_percs.groupby(['quantity', 'source', 'growth_rate_hr']):
    med = d[d['interval'] == 'median']['lower'].values[0]
    perc = d[d['interval'] == '95%']
    _df = pd.DataFrame({'quantity': g[0],
                        'source': g[1],
                        'growth_rate_hr': g[2],
                        'median_value': med,
                        '97.5%': perc['lower'].values[0],
                        '2.5%': perc['upper'].values[0]},

                       index=[0])
    ms_par_wide = pd.concat([ms_par_wide, _df], sort=False)
ms_par_wide.to_csv(
    '../../../data/mcmc/mass_spec_empirical_summaries_wide.csv', index=False)

# %%
# Summarize the empirical caluclations for the size data
post = samples.posterior['size_phiRb'].to_dataframe().reset_index()
for i in range(len(size_data)):
    post.loc[post[f'size_phiRb_dim_0'] == i,
             'source'] = size_data['source'].values[i]
    post.loc[post[f'size_phiRb_dim_0'] == i,
             'growth_rate_hr'] = size_data['growth_rate_hr'].values[i]
    post.loc[post[f'size_phiRb_dim_0'] == i,
             'surface_to_volume'] = size_data['surface_to_volume'].values[i]
    post.loc[post[f'size_phiRb_dim_0'] == i,
             'width'] = size_data['width_um'].values[i]

percs = size.viz.compute_percentiles(post, 'size_phiRb', ['source', 'width', 'surface_to_volume', 'growth_rate_hr'],
                                     lower_bounds=lowers,
                                     upper_bounds=uppers,
                                     interval_labels=labels)
percs.rename(columns={'size_phiRb': 'phiRb'}, inplace=True)

percs.to_csv(
    '../../../data/mcmc/size_data_empirical_summaries_longform.csv', index=False)

size_par_wide = pd.DataFrame([])
for g, d in percs.groupby(['quantity', 'source', 'growth_rate_hr', 'width', 'surface_to_volume']):
    med = d[d['interval'] == 'median']['lower'].values[0]
    perc = d[d['interval'] == '95%']
    _df = pd.DataFrame({'quantity': g[0],
                        'source': g[1],
                        'growth_rate_hr': g[2],
                        'width': g[3],
                        'surface_to_volume': g[4],
                        'median_value': med,
                        '97.5%': perc['lower'].values[0],
                        '2.5%': perc['upper'].values[0]}, index=[0])
    size_par_wide = pd.concat([size_par_wide, _df], sort=False)
size_par_wide.to_csv(
    '../../../data/mcmc/size_data_empirical_summaries_wide.csv', index=False)