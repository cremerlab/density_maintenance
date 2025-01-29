#%%
import numpy as np 
from collections import defaultdict
import scipy.stats
import pandas as pd 

# Load the data files
tot_rna = pd.read_csv('./raw/total_rna.csv')
tot_prot = pd.read_csv('./raw/total_protein_biuret.csv')
growth_curves = pd.read_csv('./raw/growth_curves.csv')

# Compute the growth rate from the growth curves 
rates = defaultdict(list)
for g, d in growth_curves.groupby(['carbon_source', 'replicate']):
    popt = scipy.stats.linregress(d['elapsed_time_hr'], np.log(d['od_600nm']))
    rates[g[0]].append(popt[0])

# Apply the conversion of ug_rna_per_biomass
tot_rna['ug_rna_per_biomass'] = tot_rna['od260nm'] * 31 / tot_rna['adjusted_od600nm']
tot_rna = tot_rna[['strain', 'carbon_source', 'replicate', 'ug_rna_per_biomass']]

# Load the biuret calibration file
cal = pd.read_csv('./raw/biuret_calibration_curve.csv')

# Perform a simple linear regression
popt = scipy.stats.linregress(cal['protein_conc_ug_ml'], cal['od_555nm'])
#%%
# Using the fit, compute the ug_protein_per_biomass
tot_prot['ug_protein_per_mL'] =  0.2 * (tot_prot['od555nm'] - tot_prot['od555nm_neg_control']) / (popt[0] * tot_prot['culture_volume_mL'])
tot_prot['ug_protein_per_biomass'] = tot_prot['ug_protein_per_mL'] / tot_prot['adjusted_od600nm']
tot_prot = tot_prot[['strain', 'carbon_source', 'replicate', 'ug_protein_per_biomass']]

# Merge together
merged = tot_rna.merge(tot_prot, on=['strain', 'carbon_source', 'replicate'], how='inner')

# Update growth rates. 
for k, v in rates.items():
    for i, lam in enumerate(v):
        merged.loc[(merged['carbon_source'] == k) & (merged['replicate'] == (i+1)),
                   'growth_rate_hr'] = lam
merged.to_csv('./rna_protein_per_biomass.csv', index=False)


