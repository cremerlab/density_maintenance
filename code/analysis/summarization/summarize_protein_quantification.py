# %%
import numpy as np
import pandas as pd
import scipy.stats

# Load the calibration curve data and actual measurements.
cal_data = pd.read_csv(
    '../../../data/protein_quantification/bradford_calibration_curve.csv')
prot_data = pd.read_csv(
    '../../../data/protein_quantification/wildtype_bradford_periplasmic_protein.csv')

# Perform a curve fit on the calibration data
cal_data
popt = scipy.stats.linregress(
    cal_data['protein_conc_ug_ml'], cal_data['od_595nm'])
slope, intercept = popt[:2]

# For each bradford measurement, compute the concentration of protein
prot_data['prot_ug_per_biomass'] = ((prot_data['od_595nm'] - intercept)/slope) * prot_data['dilution_factor'] * \
    prot_data['extract_volume'] / \
    (prot_data['culture_volume'] * prot_data['od_600nm'])
prot_data.to_csv(
    '../../../data/summaries/summarized_protein_measurements.csv', index=False)
