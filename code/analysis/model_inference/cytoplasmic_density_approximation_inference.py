# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()

model = cmdstanpy.CmdStanModel(
    stan_file='./cytoplasmic_density_approximation_inference.stan')

size_data = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
dna_data = pd.read_csv(
    '../../../data/literature/collated_dna_protein_ratio.csv')

# %%
phi_cyt = ms_data[ms_data['localization'] == 'cytoplasm']
phi_rib = ms_data[ms_data['localization'] == 'ribosomal sector']
N_pred = 200
pred_lam = np.linspace(0, 2.5, N_pred)
data_dict = {
    'N_size': len(size_data),
    'surface_areas': size_data['surface_area_um2'].values,
    'volume': size_data['volume_um3'].values,
    'size_lam': size_data['growth_rate_hr'].values,
    'N_dna': len(dna_data),
    'dna_protein_ratio': dna_data['DNA_protein_ratio'].values,
    'N_prot': len(prot_data),
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,
    'prot_per_cell_lam': prot_data['growth_rate_hr'].values,

    'N_ms': len(phi_cyt),
    'ms_lam': phi_cyt['growth_rate_hr'].values,
    'phi_rib': phi_rib['mass_frac'].values,
    'phi_cyt': phi_cyt['mass_frac'].values,

    'N_pred': N_pred,
    'pred_lam': pred_lam
}


_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)
# prot_data = pd.read_csv('../../data/')
# %%
post = samples.posterior['rho_cyt'].to_dataframe().reset_index()
df = pd.DataFrame([])
for g, d in post.groupby('rho_cyt_dim_0'):
    source = phi_cyt['dataset_name'].values[g]
    mean_val = np.mean(d['rho_cyt'].values)
    median_val = np.median(d['rho_cyt'].values)
    lower, upper = np.percentile(d['rho_cyt'].values, [2.5, 97.5])
    _df = pd.DataFrame({'source': source,
                        'growth_rate_hr': phi_cyt['growth_rate_hr'].values[g],
                        'mean_value': mean_val,
                        'median_value': median_val,
                        '2.5%': lower,
                        '97.5%': upper,
                        'quantity': 'rho_cyt',
                        'units': 'fg/fL'},
                       index=[0])
    df = pd.concat([df, _df], sort=False)
df.to_csv('../../../data/mcmc/approximated_cytoplasmic_density_wide.csv', index=False)

# %%

pars = ['M_prot_pred', 'DNA_pred', 'SA_pred', 'vol_pred']
df = pd.DataFrame([])
for i, p in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(post, p, f'{p}_dim_0',
                                         lower_bounds=[2.5, 12.5, 37.5, 50],
                                         upper_bounds=[97.5, 87.5, 62.5, 50],
                                         interval_labels=['95%', '75%', '25%', 'median'])
    percs['growth_rate_hr'] = [pred_lam[k] for k in percs[f'{p}_dim_0'].values]
    percs = percs[['lower', 'upper', 'quantity', 'interval', 'growth_rate_hr']]
    df = pd.concat([df, percs], sort=False)

df.to_csv('../../../data/mcmc/approximated_cytoplasmic_density_ppcs.csv', index=False)
