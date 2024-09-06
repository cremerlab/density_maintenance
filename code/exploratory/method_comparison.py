#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# Load Richa's data for comparison
richa_data = pd.read_csv('./Relative_Protein_concentrations_and_massfractions.csv')
gene_class = pd.read_csv('../../data/literature/genes_classification_all.csv')

# Define the ribosomal proteins
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

# Define the localization
locs = {'membrane': ['IM', 'LPI', 'LPO', 'OM'],
       'periplasm': ['PE']}
filt = pd.DataFrame([])
uncat = 0
for g, d in richa_data.groupby('Gene names  (primary )'):
    lcz = gene_class.loc[(gene_class['gene']==g)]
    if len(lcz) == 0:
        print(f'could not classify {g}')     
        uncat += 1
        continue
    if lcz['location'].values[0] in locs['membrane']:
        d['localization'] = 'phi_mem'
    elif lcz['location'].values[0] in locs['periplasm']:
        d['localization'] = 'phi_peri'
    elif g.lower() in ribo_prots:
        d['localization'] = 'phi_rib'
    elif g in ['lacZ']:
        d['localization'] = 'phiX-lacZ'
    elif g in ['relA']:
        d['localization'] = 'phiX-relA'
    else:
        continue 
    filt = pd.concat([filt, d])

#%%
# Set a mapper for mass fraction idx and growth rates
relative_keys = [k for k in filt.columns if ('Relative' in k) & ('log2' not in k)]
growth_rates = {'ECOR2': {'glucose': 0.8861, 'acetate': 0.2738, 'strain':'ECOR02',
                          'idx': [0,5]},
                'ECOR51': {'glucose': 0.3208, 'acetate':0.2756, 'strain':'ECOR51',
                           'idx': [1,6]},
                'ECOR63': {'glucose': 0.5048, 'acetate': 0.5526, 'strain':'ECOR63',
                           'idx': [2,7]},
                'AC': {'glucose':0.7299, 'acetate': 0.2021, 'strain':'AC1',
                       'idx': [3,8]},
                'NCM': {'glucose':0.9823, 'acetate':0.4069, 'strain':'NCM3722',
                        'idx': [4,9]}}
tidy_df = pd.DataFrame([])
for k, v in growth_rates.items():
    for i, g in enumerate(['glucose', 'acetate']):
        _df = filt[[f'Relative_concentration_{k}_{g[0].upper()+g[1:3]}', 
                    f"massfraction_{v['idx'][i]+1}", 'localization']].copy()
        _df['strain'] = v['strain']
        _df['carbon_source'] = g
        _df['growth_rate_hr'] = v[g]
        _df.rename(columns={f'Relative_concentration_{k}_{g[0].upper()+g[1:3]}': 'rel_conc',
                             f"massfraction_{v['idx'][i]+1}": 'mass_frac'}, inplace=True)     
        tidy_df = pd.concat([tidy_df, _df])

#%% 
# Sum and aggregate    
agged = tidy_df.groupby(['strain', 'carbon_source', 'localization', 'growth_rate_hr'])[['mass_frac', 'rel_conc']].sum().reset_index()

#%% 
# Load my data and lit
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']
lit_data = pd.read_csv('../../data/literature/collated_mass_fractions_empirics.csv')

#%%
fig, ax = plt.subplots(1, 3, figsize=(6,2))

for g, d in lit_data.groupby('dataset_name'):
    lit_phi_rib = d[d['localization'] == 'ribosomal sector']
    lit_phi_mem = d[d['localization'] == 'membrane']
    lit_phi_peri = d[d['localization'] == 'periplasm']
    pt = mapper[g]
    style = {'ms': 4, 'marker': pt['m'], 'markerfacecolor': pt['c'], 'markeredgecolor': 'k', 'label': g,
             'alpha':0.5}
    for i, _d in enumerate([lit_phi_rib, lit_phi_mem, lit_phi_peri]):
        ax[i].plot(_d['growth_rate_hr'], _d['mass_frac'], linestyle='none', **style)

style = {"ms": 4, "marker": 'o', 'markerfacecolor': 'w', 'markeredgecolor': cor['primary_black'],
         "markeredgewidth":1}
ax[0].plot(data['growth_rate_hr'], data['phi_rib'], 'o', label='norm (Griffin)', **style)
ax[1].plot(data['growth_rate_hr'], data['phi_mem'], 'o',**style)
ax[2].plot(data['growth_rate_hr'], data['phi_peri'], 'o',**style)


# PLot richa's data
axmapper = {'phi_rib':ax[0], 'phi_mem':ax[1], 'phi_peri':ax[2]}
style = {'ms':4, 'marker':'s', 'markerfacecolor':cor['primary_blue'], 'markeredgecolor':'k',
        'label': 'norm (Richa)', 'markeredgewidth':1}
for g, d in agged[agged['localization'].isin(axmapper.keys())].groupby('localization'):
    axmapper[g].plot(d['growth_rate_hr'], d['mass_frac'], 'o', **style)

style = {'ms':4, 'marker':'s', 'markerfacecolor':cor['primary_red'], 'markeredgecolor':'k',
        'label': 'norm (Richa)', 'markeredgewidth':1}

for g, d in agged[(agged['localization'].isin(axmapper.keys())) &
                  (agged['strain']=='NCM3722')].groupby('localization'):
    axmapper[g].plot(d['growth_rate_hr'], d['mass_frac'], 'o', **style)


for a in ax.ravel():
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\phi_{rib}$\n ribosomal mass fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$\n membrane mass fraction', fontsize=6)
ax[2].set_ylabel('$\phi_{peri}$\n periplasmic mass fraction', fontsize=6)
ax[0].set_ylim([0, 0.3])
ax[1].set_ylim([0.05, 0.2])
ax[2].set_ylim([0, 0.15])
plt.tight_layout()