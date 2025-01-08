#%%
import scipy.stats
import pandas as pd 

# Load the two data files
tot_rna = pd.read_csv('./raw/total_rna.csv')
tot_prot = pd.read_csv('./raw/total_protein_biuret.csv')

# Apply the conversion of ug_rna_per_biomass
tot_rna['ug_rna_per_biomass'] = tot_rna['od260nm'] * 31 / tot_rna['adjusted_od600nm']
tot_rna = tot_rna[['strain', 'carbon_source', 'replicate', 'ug_rna_per_biomass']]

# Load the biuret calibration file
cal = pd.read_csv('./raw/biuret_calibration_curve.csv')

# Perform a simple linear regression
popt = scipy.stats.linregress(cal['protein_conc_ug_ml'], cal['od_555nm'])

# Using the fit, compute the ug_protein_per_biomass
tot_prot['ug_protein_per_mL'] =  0.2 * (tot_prot['od555nm'] - popt[1]) / (popt[0] * tot_prot['culture_volume_mL'])
tot_prot['ug_protein_per_biomass'] = tot_prot['ug_protein_per_mL'] / tot_prot['adjusted_od600nm']
tot_prot = tot_prot[['strain', 'carbon_source', 'replicate', 'ug_protein_per_biomass']]

# Merge together
merged = tot_rna.merge(tot_prot, on=['strain', 'carbon_source', 'replicate'], how='inner')
merged.to_csv('./rna_protein_per_biomass.csv', index=False)


