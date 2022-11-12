# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz

# Load the data
data = pd.read_csv(
    '../../data/protein_quantification/wildtype_protein_quantification.csv')
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/simple_prot_quantification.stan')

percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']

samples_df, summary_df = [], []
for g, d in data.groupby(['carbon_source']):
    data_dict = {'N': len(d),
                 'm_cyto': d['m_cyto'].values.astype(float),
                 'm_peri': d['m_peri'].values.astype(float)}
    samples = model.sample(data=data_dict)
    samples = arviz.from_cmdstanpy(samples)
    samp_df = samples.posterior.to_dataframe().reset_index()
    samp_df['phi_peri'] = samp_df['m_peri_mu'].values / \
        (samp_df['m_peri_mu'] + samp_df['m_cyto_mu'])
    for var in ['m_cyto_mu', 'm_peri_mu', 'phi_peri']:
        _d = samp_df[var].values
        _percs = np.percentile(_d, percs)
        _df = pd.DataFrame([_percs], columns=perc_cols)
        _df['median'] = np.median(_d)
        _df['carbon_source'] = g
        _df['parameter'] = var.split('_mu')[0]
        summary_df.append(_df)
    samples_df.append(samp_df)

summary_df = pd.concat(summary_df, sort=False)
samples_df = pd.concat(samples_df, sort=False)
summary_df.to_csv(
    '../../data/protein_quantification/protein_quantification_summary.csv', index=False)
samples_df.to_csv(
    '../../data/protein_quantification/protein_quantification_samples.csv', index=False)
