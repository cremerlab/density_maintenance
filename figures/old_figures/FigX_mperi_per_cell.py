# %%
import scipy.stats
import numpy as np
import pandas as pd
import size.viz
cor, pal = size.viz.matplotlib_style()

growth_data = pd.read_csv(
    '../data/summaries/summarized_growth_measurements.csv')
size_data = pd.read_csv(
    '../data/summaries/summarized_size_measurements.csv')
prot_data = pd.read_csv(
    '../data/summaries/summarized_protein_measurements.csv')
params = pd.read_csv('../data/mcmc/parameter_percentiles.csv')

# Group the datasets
mass_spec = pd.read_csv('../data/literature/compiled_mass_fractions.csv')
mass_spec = mass_spec[mass_spec['periplasmic'] == True]
# mass_spec = mass_spec.groupby(['dataset_name', 'condition', 'growth_rate_hr'])[
# 'mass_frac'].agg(('mean', 'sem')).reset_index()

# Group the measurements
growth_data = growth_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL'])[
    'growth_rate_hr'].agg(('mean', 'sem'))
size_data = size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc'])[
    'growth_rate_hr'].agg(('mean', 'sem'))
prot_data = prot_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_ml'])[
    'growth_rate_hr'].agg(('mean', 'sem'))
flow_data = flow_data.groupby(['strain', 'carbon_source']).agg(('mean', 'sem'))

# Link growth rates to various data

# for g, d in size_data[(size_data['strain'] == 'wildtype') &
#   (size_data['overexpression'] == 'none') &
#   (size_data['inducer_conc'] == 0)].groupby(['carbon_source']):


# flow_fit = scipy.stats.linregress
