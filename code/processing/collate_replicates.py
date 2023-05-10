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
growth_data.rename(columns={'temperature_C': 'temperature'}, inplace=True)
size_data.rename(columns={'temperature_C': 'temperature'}, inplace=True)
reps = pd.read_csv('../../data/replicate_analysis.csv')
flow = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
flow = flow[(flow['strain'] == 'wildtype') & (flow['carbon_source'] != 'LB')]
flow = flow.groupby(['date', 'run_no', 'carbon_source'])[
    'cells_per_biomass'].mean().reset_index()
biuret = pd.read_csv(
    '../../data/protein_quantification/biuret_total_protein.csv')

biuret = biuret[biuret['growth_rate_hr'] < 1.5]

# Use only the accepted perturbations
reps = reps[reps['accepted'] == True]
reps['idx'] = np.arange(len(reps)) + 1

carbon_label = {'acetate': 1, 'sorbitol': 2,
                'glycerol': 3, 'glucose': 4, 'glucoseCAA': 5, 'LB': 6}
reps['carbon_idx'] = [carbon_label[g] for g in reps['carbon_source'].values]

# Define the carbon idxs
wt_carbon_label = {}
for g, d in reps[(reps['strain'] == 'wildtype') & (reps['overexpression'] == 'none') &
                 (reps['inducer_conc'] == 0) & (reps['temperature'] == 37)].groupby(['carbon_source', 'idx']):
    wt_carbon_label[g[0]] = g[1]
# %%
# Ensure there's no listed inducer witha a concentration 0
for d in [size_data, brad_data, growth_data]:
    d.loc[(d['inducer'] != 'none') & (
        d['inducer_conc'] == 0), 'inducer'] = 'none'

valid_dfs = []

keys = ['strain', 'carbon_source', 'overexpression',
        'inducer', 'inducer_conc', 'temperature']
for i, d in enumerate([size_data, growth_data, brad_data]):
    _df = pd.DataFrame([])
    for j in range(len(reps)):
        _s, _c, _o, _ind, _indc, _tem = reps.iloc[j][keys].values

        _d = d[(d['strain'] == _s) & (d['carbon_source'] == _c) &
               (d['overexpression'] == _o) & (d['inducer'] == _ind) &
               (d['inducer_conc'] == _indc) & (d['temperature'] == _tem)].copy()
        if len(_d) > 0:
            idx = j + 1  # len(reps) - (j + 1)
            _d['perturbation_idx'] = idx
            wt_idx = reps.iloc[j]['carbon_source']
            _d['wt_carbon_idx'] = int(
                wt_carbon_label[reps.iloc[j]['carbon_source']])
            _d['carbon_idx'] = int(
                carbon_label[reps.iloc[j]['carbon_source']])

            _df = pd.concat([_df, _d], sort=False)
        else:
            print(reps.iloc[j][keys].values)
    valid_dfs.append(_df)

for d in [flow, biuret]:
    for k, v in wt_carbon_label.items():
        d.loc[d['carbon_source'] == k, 'wt_carbon_idx'] = int(v)
    d.dropna(inplace=True)
    d['cond_idx'] = d.groupby(['wt_carbon_idx']).ngroup() + 1
    d['carbon_idx'] = [carbon_label[g] for g in d['carbon_source'].values]
    d['wt_carbon_idx'] = d['wt_carbon_idx'].values.astype(int)
    valid_dfs.append(d)

# Add specific mappers to the brad data
valid_dfs[2]['idx'] = valid_dfs[2].groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']).ngroup() + 1


# %%
# Save the pruned and labeled dataframes to disk.
names = ['size', 'growth', 'protein', 'flow', 'total_protein']
for i, nom in enumerate(names):
    valid_dfs[i].to_csv(
        f'../../data/processed/labeled_{nom}_measurements.csv', index=False)
