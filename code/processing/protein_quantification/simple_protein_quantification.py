# %%
import numpy as np
import pandas as pd
import scipy.stats

# Load the calibration curve and the raw measurements
calib = pd.read_csv(
    '../../../data/protein_quantification/bradford_calibration_curve.csv')
data = pd.read_csv(
    '../../../data/protein_quantification/bradford_periplasmic_protein.csv')

# Do a simple linear regression on the pooled calibration curve data
popt = scipy.stats.linregress(calib['protein_conc_ug_ml'], calib['od_595nm'])
slope = popt[0]
intercept = popt[1]
# %%
# Convert the measured OD595 in the data to ug protein
data['ug_protein'] = (data['od_595nm'] - intercept) / slope
data['m_peri_per_biomass'] = data['ug_protein'] * data['extraction_volume_ml'] * \
    data['dilution_factor'] / \
    (data['od_600nm'] * data['culture_volume_ml'])
# simple_quantification = data[[
# 'strain', 'carbon_source', 'date', 'm_peri_per_biomass']]
data.to_csv(
    '../../../data/protein_quantification/periplasmic_protein_simple_quantification.csv', index=False)
