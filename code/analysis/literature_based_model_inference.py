# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
import scipy.stats
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
total = ms_data[(ms_data['localization'] == 'cytoplasm') | (
    ms_data['localization'] == 'envelope')].groupby(
        ['dataset_name', 'condition', 'growth_rate_hr', 'volume']
)['mass_fg'].sum().reset_index()

#
vol_scaling_pars = [['w_min', 'w_slope'], ['tau', 'C_0'],]
corner_pars = [['alpha', 'rho_prot_min', 'rho_prot_slope',
               'm_peri_mu', 'w_sigma', 'ell_sigma', 'vol_sigma',
                'phi_peri_sigma', 'rho_mem_mu'],
               ['alpha', 'rho_prot_min', 'rho_prot_slope',
               'm_peri_mu', 'w_sigma', 'ell_sigma', 'vol_sigma',
                'phi_peri_sigma', 'phi_mem_sigma', 'phi_mem_mu']]

size_pars = ['w', 'ell', 'vol']
model_pars = [['phi_peri', 'phi_mem', 'rho_prot',
              'rho_peri', 'rel_phi',
               'rho_mem', 'alpha', 'm_peri_mu', 'phi_mem_mu'],
              ['phi_peri', 'phi_mem', 'rho_prot', 'alpha', 'phi_mem',
              'rho_peri', 'm_peri', 'rel_phi', 'rho_mem']]
model_kde_pars = [['rho_prot', 'rho_mem_mu', 'alpha', 'm_peri_mu'],
                  ['rho_prot', 'phi_mem_mu', 'alpha', 'm_peri_mu']]


prot_pars = ['prot_per_cell']
model_desc = ['const_rho_mem', 'const_phi_mem']
vol_desc = ['linear_width', 'smk']

# Compute the percentiles for the parameters
upper_percs = [97.5, 87.5, 75, 62.5, 55, 50]
lower_percs = [2.5, 12.5, 25, 37.5, 45, 50]
labels = ['95%', '75%', '50%', '25%', '10%', 'median']
kwargs = {'lower_bounds': lower_percs,
          'upper_bounds': upper_percs, 'interval_labels': labels}

perc_df = pd.DataFrame([])
kde_df = pd.DataFrame([])
for i in range(2):
    for j in range(2):
        N_sim = 1000
        lam_sim = np.linspace(0, 3.1, N_sim)
        width_sim = np.linspace(0.45, 1.5, N_sim)
        # Assemble the data dictionary
        data_dict = {
            'smk': i,
            'const_phi_mem': j,
            'N_size': len(size_data),
            'N_prot': len(prot_data),
            'N_mass_spec': len(membrane),
            'N_sim': N_sim,
            'lam_sim': lam_sim,
            'width_sim': width_sim,
            'delta': 0.0249,
            'widths': size_data['width_um'].values.astype(float),
            'lengths': size_data['length_um'].values.astype(float),
            'volumes': size_data['volume_um3'].values.astype(float),
            'size_lam': size_data['growth_rate_hr'].values.astype(float),
            'prot_per_cell': prot_data['fg_protein_per_cell'].values.astype(float),
            'prot_lam': prot_data['growth_rate_hr'].values.astype(float),
            'phi_mem': membrane['mass_frac'].values.astype(float),
            'rho_mem_meas': membrane['mass_fg'].values.astype(float) / (2 * membrane['surface_to_volume'].values.astype(float) * membrane['volume'].astype(float)),
            'rho_prot_meas': prot_data['density'].values.astype(float),
            'rho_prot_meas_ms': total['mass_fg'].values.astype(float) / total['volume'].values.astype(float),
            'phi_peri': periplasm['mass_frac'].values.astype(float),
            'm_peri_meas': periplasm['mass_fg'].values.astype(float),
            'ms_lam': periplasm['growth_rate_hr'].values.astype(float)}

        # Sample the model
        _samples = model.sample(data=data_dict)
        samples = az.from_cmdstanpy(_samples)

        # Finished sampling
        fig = plt.figure(figsize=(10, 10))
        fig = corner(samples, group='posterior', var_names=vol_scaling_pars[i] + corner_pars[j], fig=fig,
                     hist_kwargs={'lw': 1}, plot_contours=False, plot_density=False, data_kwargs={'ms': 1},
                     divergences=True, divergences_kwargs={'color': cor['primary_red'], 'ms': 1, 'markeredgewidth': 0})
        for a in fig.axes:
            a.grid(False)
        plt.savefig(
            f'../../figures/mcmc/model_diagnostics/model{j+1}_smk{i}_corner.pdf', bbox_inches='tight')
        plt.close()

        post = samples.posterior[vol_scaling_pars[i] +
                                 corner_pars[j]].to_dataframe().reset_index()
        post['idx'] = 1
        percs = size.viz.compute_percentiles(
            post, vol_scaling_pars[i] + corner_pars[j], 'idx', **kwargs)
        percs.to_csv(
            f'../../data/mcmc/literature_model{j+1}_smk{i}_parameter_percs.csv', index=False)

        # Compute the percentiles ppc
        datasets = size_data['source'].values
        growth_rates = size_data['growth_rate_hr'].values
        _perc_df = pd.DataFrame([])
        for _j in range(2):
            print("Processing...")
            if _j == 0:
                suff = 'rep'
            else:
                suff = 'sim'
            for k, pars in enumerate([size_pars, model_pars[j], prot_pars]):
                for p in pars:
                    print("Computing ppcs...")
                    if k == 1:
                        p = p.split('_mu')[0]
                    ppc = samples.posterior[f'{p}_{suff}'].to_dataframe(
                    ).reset_index()
                    ppc.dropna(inplace=True)
                    percs = size.viz.compute_percentiles(
                        ppc, f'{p}_{suff}', f'{p}_{suff}_dim_0', **kwargs)
                    for ell, (lam, width) in enumerate(zip(lam_sim, width_sim)):
                        percs.loc[percs[f'{p}_{suff}_dim_0']
                                  == ell, 'growth_rate_hr'] = lam
                        percs.loc[percs[f'{p}_{suff}_dim_0']
                                  == ell, 'width'] = width
                    percs.drop(columns=[f'{p}_{suff}_dim_0'], inplace=True)
                    _perc_df = pd.concat([_perc_df, percs], sort=False)
        _perc_df['volume_scale'] = vol_desc[i]
        _perc_df['model'] = model_desc[j]
        perc_df = pd.concat([perc_df, _perc_df], sort=False)

        # Compute a kernel density estimate over all of the model pars.
        for _, p in enumerate(model_kde_pars[j]):
            post = samples.posterior[p].to_dataframe().reset_index()
            prange = np.linspace(
                0.75 * post[p].min(), 1.25 * post[p].max(), 300)
            kernel = scipy.stats.gaussian_kde(post[p].values)
            kde = kernel(prange)
            kde *= kde.sum()**-1
            _df = pd.DataFrame(
                np.array([prange, kde]).T, columns=['value', 'kde'])
            _df['parameter'] = p
            _df['model'] = model_desc[j]
            _df['volume_scale'] = vol_desc[i]
            kde_df = pd.concat([kde_df, _df], sort=False)

perc_df.to_csv(
    f'../../data/mcmc/literature_model_params.csv', index=False)
kde_df.to_csv('../../data/mcmc/literature_model_params_kde.csv', index=False)
print("Done!")
