# %%

import numpy as np
import pandas as pd
data = pd.read_csv(
    '../../../data/mass_spectrometry/periplasm_fractionation_confirmation_raw.csv')
data.dropna(inplace=True)
data['cyto_enrichment'] = data['cyto_intensity'] / data['proteome_intensity']
data['peri_enrichment'] = data['peri_intensity'] / data['proteome_intensity']
data['relative_signal'] = data['peri_intensity'] / data['cyto_intensity']
data['fractional_peri'] = data['peri_intensity'] / data['proteome_intensity'].sum()
data['fractional_cyto'] = data['cyto_intensity'] / data['proteome_intensity'].sum()
data = data[['gene', 'localization', 'cyto_enrichment',
             'peri_enrichment', 'relative_signal', 'fractional_peri', 'fractional_cyto']]
data.dropna(inplace=True)
data.to_csv(
    '../../../data/mass_spectrometry/periplasm_fractionation_confirmation.csv', index=False)