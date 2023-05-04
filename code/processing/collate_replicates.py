# %%
import numpy as np
import pandas as pd

# Load the various datasets
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
brad_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein.csv')
growth_data = pd.read_csv(
    '../../data/summaries/summarized_growth_measurements.csv')
reps = pd.read_csv('../../data/replicate_analysis.csv')
flow = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
flow = flow[flow['strain'] == 'wildtype']
biuret = pd.read_csv(
    '../../data/protein_quantification/biuret_total_protein.csv')

# Define the carbon idxs
carbon_label = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3, 'glucose': 4,
                'glucoseCAA': 5, 'LB': 21}

# Ensure there's no listed inducer witha a concentration 0
for d in [size_data, brad_data, growth_data]:
    d.loc[(d['inducer'] != 'none') & (
        d['inducer_conc'] == 0), 'inducer'] = 'none'

# Use only the accepted perturbations
reps = reps[reps['accepted'] == True]

valid_dfs = []

keys = ['strain', 'carbon_source', 'overexpression', 'inducer_conc']
for i, d in enumerate([size_data, growth_data, brad_data]):
    _df = pd.DataFrame([])
    for j in range(len(reps)):
        _s, _c, _o, _ind = reps.iloc[j][keys].values
        _d = d[(d['strain'] == _s) & (d['carbon_source'] == _c) &
               (d['overexpression'] == _o) & (d['inducer_conc'] == _ind)].copy()
        if len(_d) > 0:
            idx = len(reps) - (j + 1)
            if idx == 0:
                idx = len(reps)
            _d['perturbation_idx'] = idx
            _d['wt_carbon_idx'] = int(
                carbon_label[reps.iloc[j]['carbon_source']])
            _df = pd.concat([_df, _d], sort=False)
    valid_dfs.append(_df)

for d in [flow, biuret]:
    for k, v in carbon_label.items():
        d.loc[d['carbon_source'] == k, 'wt_carbon_idx'] = int(v)
    d.dropna(inplace=True)
    d['wt_carbon_idx'] = d['wt_carbon_idx'].values.astype(int)
    valid_dfs.append(d)

# Save the pruned and labeled dataframes to disk.
names = ['size', 'growth', 'protein', 'flow', 'total_protein']
for i, nom in enumerate(names):
    valid_dfs[i].to_csv(
        f'../../data/processed/labeled_{nom}_measurements.csv', index=False)
