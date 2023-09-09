# %%
import numpy as np
import pandas as pd

prot = pd.read_csv('../../../data/summaries/summarized_total_protein.csv')
rna = pd.read_csv('../../../data/summaries/summarized_total_rna.csv')

collated = pd.DataFrame([])

for g, d in prot.groupby(['strain', 'overexpression', 'carbon_source', 'inducer_conc', 'replicate']):
    _rna = rna[(rna['strain'] == g[0]) &
               (rna['overexpression'] == g[1]) &
               (rna['carbon_source'] == g[2]) &
               (rna['inducer_conc'] == g[3]) &
               (rna['replicate'] == g[4])]
    d = d.copy()

    d['ug_rna_per_biomass'] = _rna['ug_rna_per_biomass'].values
    d['phiRb'] = 0.4558 * d['ug_rna_per_biomass'].values[0] / \
        d['ug_prot_per_biomass'].values[0]
    d = d[['strain', 'overexpression', 'carbon_source', 'inducer_conc', 'replicate',
           'ug_rna_per_biomass', 'ug_prot_per_biomass', 'phiRb']]
    collated = pd.concat([collated, d])

collated.to_csv('../../../data/summaries/merged_rna_protein_measurements.csv')
