# %%
import numpy as np
import pandas as pd

rna = pd.read_csv('../../../data/protein_quantification/total_rna.csv')
rna['ug_rna_per_biomass'] = rna['od260nm'] * 31 / rna['adjusted_od600nm']
rna = rna[rna['valid'] == True]
rna = rna[['strain', 'overexpression', 'inducer_conc',
           'carbon_source', 'replicate', 'ug_rna_per_biomass']]
rna.to_csv('../../../data/summaries/summarized_total_rna.csv', index=False)
