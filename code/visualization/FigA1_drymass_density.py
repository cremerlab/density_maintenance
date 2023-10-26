# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import size.viz
cor, pal = size.viz.matplotlib_style()
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
drymass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')
prot = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
dna = pd.read_csv('../../data/literature/collated_dna_protein_ratio.csv')
vol_data = pd.read_csv(
    '../../data/literature/collated_literature_volume_data.csv')

dna_const = dna['DNA_protein_ratio'].values.mean()

# Do an empirical exponential fit to the protein data
prot_popt = scipy.stats.linregress(
    prot['growth_rate_hr'], np.log(prot['fg_protein_per_cell']))
vol_popt = scipy.stats.linregress(
    vol_data['growth_rate_hr'], np.log(vol_data['volume_um3']))

est_ms_rho_cyt = pd.DataFrame([])
for g, d in ms_data.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    cyt = d[d['localization'] == 'cytoplasm']
    rib = d[d['localization'] == 'ribosomal sector']
    M_prot = np.exp(prot_popt[1] + prot_popt[0] *
                    cyt['growth_rate_hr'].values[0])
    V = np.exp(vol_popt[1] + vol_popt[0] * cyt['growth_rate_hr'].values[0])
    M_RNA = rib['mass_frac'].values[0] * M_prot / 0.4558
    M_DNA = M_prot * dna_const
    tot_drymass = cyt['mass_frac'].values[0] * M_prot + M_RNA + M_DNA
    density = tot_drymass / V
    _df = pd.DataFrame({'source': g[0],
                        'growth_rate_hr': g[1],
                        'condition': g[2],
                        'total_protein': M_prot,
                        'cytoplasmic_protein': cyt['mass_frac'] * M_prot,
                        'total_RNA': M_RNA,
                        'total_DNA': M_DNA,
                        'total_drymass': tot_drymass,
                        'rho_cyt': density}, index=[0])
    est_ms_rho_cyt = pd.concat([est_ms_rho_cyt, _df], sort=False)

# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_ylim(100, 450)
for g, d in drymass.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **fmt)

for g, d in est_ms_rho_cyt.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['primary_red']
    ax.plot(d['growth_rate_hr'], d['rho_cyt'], **fmt)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('mass density [fg / fL]', fontsize=6)
# ax.legend()
plt.savefig('../../figures/FigA1_cytoplasmic_drymass_density.pdf',
            bbox_inches='tight')
