# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
import arviz as az
import cmdstanpy
import size.viz
import tqdm
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Compile the model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/literature_based_model_inference.stan')

# Load the data sets
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
prot_data = pd.read_csv(
    '../../data/literature/collated_total_protein_density.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')

# Aggregate the mass spec data
membrane = ms_data[ms_data['localization'] == 'membrane']
periplasm = ms_data[ms_data['localization'] == 'periplasm']

#
corner_pars = [['w_min', 'w_slope', 'alpha', 'rho_prot_min', 'rho_prot_slope', 'rho_prot_sigma',
               'm_peri_mu', 'w_sigma', 'ell_sigma', 'vol_sigma',
                'phi_peri_sigma', 'rho_mem_mu'],
               ['w_min', 'w_slope', 'alpha', 'rho_prot_min', 'rho_prot_slope', 'rho_prot_sigma',
               'm_peri_mu', 'w_sigma', 'ell_sigma', 'vol_sigma',
                'phi_peri_sigma', 'phi_mem_sigma', 'phi_mem_mu']]

size_pars = ['w', 'ell', 'vol']
model_pars = [['phi_peri', 'phi_mem', 'rho_prot',
              'rho_peri', 'm_peri',  'rel_phi',
               'rho_mem'],
              ['phi_peri', 'phi_mem', 'rho_prot',
              'rho_peri', 'm_peri', 'rel_phi', 'phi_mem',
               'rho_mem']]


prot_pars = ['prot_per_cell']
model_desc = ['const_rho_mem', 'const_phi_mem']
# Compute the percentiles for the parameters
upper_percs = [97.5, 87.5, 75, 62.5, 55, 50]
lower_percs = [2.5, 12.5, 25, 37.5, 45, 50]
labels = ['95%', '75%', '50%', '25%', '10%', 'median']
kwargs = {'lower_bounds': lower_percs,
          'upper_bounds': upper_percs, 'interval_labels': labels}

perc_df = pd.DataFrame([])
for i in tqdm.tqdm(range(2)):
    N_sim = 1000
    lam_sim = np.linspace(0, 3.1, N_sim)
    # Assemble the data dictionary
    data_dict = {
        'N_size': len(size_data),
        'N_prot': len(prot_data),
        'N_mass_spec': len(membrane),
        'N_sim': N_sim,
        'lam_sim': lam_sim,
        'delta': 0.0249,
        'const_phi_mem': i,
        'widths': size_data['width_um'].values.astype(float),
        'lengths': size_data['length_um'].values.astype(float),
        'volumes': size_data['volume_um3'].values.astype(float),
        'size_lam': size_data['growth_rate_hr'].values.astype(float),
        'prot_per_cell': prot_data['fg_protein_per_cell'].values.astype(float),
        'prot_lam': prot_data['growth_rate_hr'].values.astype(float),
        'phi_mem': membrane['mass_frac'].values.astype(float),
        'rho_mem_meas': membrane['mass_fg'].values.astype(float) / (membrane['surface_to_volume'].values.astype(float) * membrane['volume'].astype(float)),
        'rho_prot_meas': prot_data['density'].values.astype(float),
        'phi_peri': periplasm['mass_frac'].values.astype(float),
        'm_peri_meas': periplasm['mass_fg'].values.astype(float),
        'ms_lam': periplasm['growth_rate_hr'].values.astype(float)}

    # Sample the model
    _samples = model.sample(data=data_dict)  # , show_console=True)
    samples = az.from_cmdstanpy(_samples)
    # Finished sampling
    fig = plt.figure(figsize=(10, 10))
    fig = corner(samples, group='posterior', var_names=corner_pars[i], fig=fig,
                 hist_kwargs={'lw': 1}, plot_contours=False, plot_density=False, data_kwargs={'ms': 1},
                 divergences=True, divergences_kwargs={'color': cor['primary_red'], 'ms': 1, 'markeredgewidth': 0})
    for a in fig.axes:
        a.grid(False)
    plt.savefig(
        f'../../figures/mcmc/model_diagnostics/model{i+1}_corner.pdf', bbox_inches='tight')
    plt.close()

    post = samples.posterior[corner_pars[i]].to_dataframe().reset_index()
    post['idx'] = 1
    percs = size.viz.compute_percentiles(post, corner_pars[i], 'idx', **kwargs)
    percs.to_csv(
        f'../../data/mcmc/literature_model{i+1}_parameter_percs.csv', index=False)

    # Compute the percentiles for the size ppc
    datasets = size_data['source'].values
    growth_rates = size_data['growth_rate_hr'].values
    _perc_df = pd.DataFrame([])
    for j in range(2):
        print("Processing...")
        if j == 0:
            suff = 'rep'
        else:
            suff = 'sim'
        for k, pars in enumerate([size_pars, model_pars[i], prot_pars]):
            for p in pars:
                print("Computing ppcs...")
                p = p.split('_mu')[0]
                ppc = samples.posterior[f'{p}_{suff}'].to_dataframe(
                ).reset_index()
                percs = size.viz.compute_percentiles(
                    ppc, f'{p}_{suff}', f'{p}_{suff}_dim_0', **kwargs)
                for ell, lam in enumerate(lam_sim):
                    percs.loc[percs[f'{p}_{suff}_dim_0']
                              == ell, 'growth_rate_hr'] = lam
                percs.drop(columns=[f'{p}_{suff}_dim_0'], inplace=True)
                _perc_df = pd.concat([_perc_df, percs], sort=False)
    _perc_df['model'] = model_desc[i]
    perc_df = pd.concat([perc_df, _perc_df], sort=False)
perc_df.to_csv(
    f'../../data/mcmc/literature_model_params.csv', index=False)
print("Done!")
