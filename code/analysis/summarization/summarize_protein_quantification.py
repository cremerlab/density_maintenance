# %%
import numpy as np
import pandas as pd
import scipy.stats

# Load the calibration curve data and actual measurements.
cal_data = pd.read_csv(
    '../../../data/protein_quantification/bradford_calibration_curve.csv')
prot_data = pd.read_csv(
    '../../../data/protein_quantification/bradford_periplasmic_protein.csv')
prot_data = prot_data[prot_data['strain'].isin(
    ['lpp14', 'malE-rbsB-fliC-KO', 'wildtype'])]

# Load the growth rate data and the literature measurements to do a simple regression
# to figure total protein as a function of growth rate
growth_rates = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')
growth_rates = growth_rates.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc'])[
    'growth_rate_hr'].mean().reset_index()
lit_data = pd.read_csv(
    '../../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')
# %%
# Perform a fit on the literature data, excluding RDM, for protein mass per od
lit_data = lit_data[lit_data['medium'] != 'RDM']
popt = scipy.stats.linregress(
    lit_data['growth_rate_hr'], np.log(lit_data['protein_mass_ug']))

# Given the fit, compute the total protein per OD given growth rates

growth_rates['protein_mass_ug'] = np.exp(
    popt[1] + popt[0] * growth_rates['growth_rate_hr'])

# %%
plt.plot(lit_data['growth_rate_hr'], lit_data['protein_mass_ug'], 'o')
lam = np.linspace(0.2, 1.5, 100)
fit = np.exp(popt[1] + popt[0] * lam)
plt.semilogy(lam, fit, 'k-')

#

# %%
# Perform a curve fit on the calibration data
cal_data
popt = scipy.stats.linregress(
    cal_data['protein_conc_ug_ml'], cal_data['od_595nm'])
slope, intercept = popt[:2]

# For each bradford measurement, compute the concentration of protein
prot_data['prot_ug_per_biomass'] = ((prot_data['od_595nm'] - intercept)/slope) * prot_data['dilution_factor'] * \
    prot_data['extraction_volume_mL'] / \
    (prot_data['culture_volume_mL'] * prot_data['od_600nm'])

dfs = []
for g, d in prot_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']):
    d = d.copy()
    lam = growth_rates[(growth_rates['strain'] == g[0]) &
                       (growth_rates['carbon_source'] == g[1]) &
                       (growth_rates['overexpression'] == g[2]) &
                       (growth_rates['inducer_conc'] == g[3])]
    if len(lam) != 0:
        print(g)
        d['tot_prot_ug_per_biomass'] = lam['protein_mass_ug'].values[0]
    dfs.append(d)
prot_data = pd.concat(dfs, sort=False)
prot_data['mass_frac'] = prot_data['prot_ug_per_biomass'] / \
    prot_data['tot_prot_ug_per_biomass']
prot_data
prot_data.to_csv(
    '../../../data/summaries/summarized_protein_measurements.csv', index=False)
