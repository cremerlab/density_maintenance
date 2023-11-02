# %%
import numpy as np
import pandas as pd
import scipy.stats
brad_cal = pd.read_csv(
    '../../../data/protein_quantification/bradford_calibration_curve.csv')
biuret_cal = pd.read_csv(
    '../../../data/protein_quantification/biuret_calibration_curve.csv')
bca_cal = pd.read_csv(
    '../../../data/protein_quantification/bca_calibration_curve.csv')

# %%
# Compute the calibration curves
brad_popt = scipy.stats.linregress(
    brad_cal['protein_conc_ug_ml'], brad_cal['od_595nm'])
biuret_popt = scipy.stats.linregress(
    biuret_cal['protein_conc_ug_ml'], biuret_cal['od_555nm'])
bca_popt = scipy.stats.linregress(
    bca_cal['protein_conc_ug_ml'], bca_cal['od_562nm'])
popts = {'bradford':brad_popt, 'biuret':biuret_popt, 'BCA':bca_popt}
comp_quant = pd.DataFrame([])
for g, d in comparison.groupby('method'):
    popt = popts[g]
    d['conc'] = d['dilution_factor'] * (d['od'])/(popt[0])
    if g == 'BCA':
        _d = d[d['dilution_factor'] <= 160]
        _d['dilution_factor'] = _d['dilution_factor'] / 10 
    else:
        _d = d[d['dilution_factor'] <= 16]
    comp_quant = pd.concat([comp_quant, _d])
# %%
# Periplasmic protein processing
peri_data = pd.read_csv(
    '../../../data/protein_quantification/bradford_periplasmic_protein.csv')
peri_data['conv_factor'] = peri_data['od_600nm'] * peri_data['dilution_factor'] * \
    peri_data['extraction_volume_mL'] / peri_data['culture_volume_mL']
peri_data['protein_conc_ug_mL'] = (
    peri_data['od_595nm'] - brad_popt[1]) / brad_popt[0]
peri_data['ug_prot_per_biomass'] = np.round(
    peri_data['protein_conc_ug_mL'] / peri_data['conv_factor'], decimals=3)
peri_data['inducer_conc'] = peri_data['inducer_conc'].values.astype(float)
peri_data['temperature'] = peri_data['temperature'].values.astype(float)
peri_data = peri_data[peri_data['strain'].isin(['wildtype', 'lpp14']) &
                      peri_data['overexpression'].isin(['none', 'relA', 'meshI']) &
                      peri_data['inducer_conc'].isin([0, 1, 100]) &
                      (peri_data['temperature'] == 37.0)]
peri_data = peri_data[['strain', 'carbon_source', 'overexpression',
                       'temperature', 'inducer_conc', 'ug_prot_per_biomass']]

peri_data.to_csv(
    '../../../data/summaries/summarized_periplasmic_protein.csv', index=False)


# Total protein processing
tot_prot = pd.read_csv(
    '../../../data/protein_quantification/total_protein_biuret.csv')
tot_prot['protein_conc_ug_ml'] = 0.2 * (
    tot_prot['adjusted_od555nm'] - biuret_popt[1]) / (biuret_popt[0] * tot_prot['culture_volume_mL'])
tot_prot['ug_prot_per_biomass'] = tot_prot['protein_conc_ug_ml'].values / \
    tot_prot['adjusted_od600nm']
tot_prot = tot_prot[tot_prot['valid'] == True]
tot_prot = tot_prot[['strain', 'overexpression', 'replicate', 'inducer_conc', 'carbon_source',
                     'ug_prot_per_biomass']]
tot_prot.to_csv(
    '../../../data/summaries/summarized_total_protein.csv', index=False)

# Membrane protein processing
mem_prot = pd.read_csv(
    '../../../data/protein_quantification/collated_BCA_measurements.csv')
mem_prot = mem_prot[mem_prot['fraction'] == 'membrane']
mem_prot['conv_factor'] = (mem_prot['od600nm'].values * mem_prot['culture_volume_mL'] /
                           (mem_prot['dilution_factor'].values * mem_prot['extraction_volume']))
# mem_prot['protein_conc_ug_ml'] = mem_prot['']
mem_prot['protein_conc_ug_ml'] = (
    mem_prot['od562nm'] - bca_popt[1]) / (bca_popt[0])
mem_prot['ug_prot_per_biomass'] = mem_prot['protein_conc_ug_ml'] / \
    mem_prot['conv_factor']
mem_prot = mem_prot[['strain', 'carbon_source', 'overexpression',
                     'inducer_conc', 'biological_replicate', 'ug_prot_per_biomass']]
mem_prot = mem_prot.groupby(['strain', 'carbon_source', 'overexpression',
                            'inducer_conc', 'biological_replicate']).mean().reset_index()
mem_prot.dropna(inplace=True)
mem_prot.to_csv(
    '../../../data/summaries/summarized_membrane_protein.csv', index=False)
