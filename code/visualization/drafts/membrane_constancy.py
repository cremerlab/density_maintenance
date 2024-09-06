# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
cor, pal = size.viz.matplotlib_style()
# mapper = size.viz.load_markercolors()
mass_spec = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
mass_spec_inf = pd.read_csv('../../data/mcmc/mass_spec_percentiles.csv')
mem_prots = mass_spec[mass_spec['go_terms'].str.contains('GO:0005886')]
mem_prots = mem_prots.groupby(['dataset_name', 'growth_rate_hr',
                               'condition'])['mass_frac'].sum().reset_index()
peri_prots = mass_spec[mass_spec['periplasm'] == True].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()

markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h',
           'p', 'P', '*', 'o', '8', 'd', '>', 'v', '<', '^']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(mem_prots['dataset_name'].unique())
# for n in mem_prots['dataset_name'].unique():
# names.append(n)
mapper = {n: {'m': m, 'c': c} for n, m, c in zip(names, markers, cors)}

# Assemble a dataframe containing membrane protein density
areal_density_df = pd.DataFrame([])

for g, d in mem_prots.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    sa = mass_spec_inf[(mass_spec_inf['dataset_name'] == g[0]) & (mass_spec_inf['growth_rate_hr'] == g[1]) &
                       (mass_spec_inf['condition'] == g[2]) &
                       (mass_spec_inf['interval'] == 'median') & (mass_spec_inf['quantity'] == 'mass_spec_sa')]['lower'].values[0]
    total_prot = mass_spec_inf[(mass_spec_inf['dataset_name'] == g[0]) & (mass_spec_inf['growth_rate_hr'] == g[1]) &
                               (mass_spec_inf['condition'] == g[2]) &
                               (mass_spec_inf['interval'] == 'median') & (mass_spec_inf['quantity'] == 'mass_spec_tot_prot_per_cell')]['lower'].values[0]
    rho = d['mass_frac'].values[0] * total_prot / sa
    prot_mass = d['mass_frac'].values[0] * total_prot
    _df = pd.DataFrame(np.array([rho]).T, columns=['areal_density'])
    _df['prot_mass'] = prot_mass
    _df['dataset_name'] = g[0]
    _df['condition'] = g[2]
    _df['growth_rate_hr'] = g[1]
    _df['total'] = total_prot
    areal_density_df = pd.concat([areal_density_df, _df], sort=False)


# %%

# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
for g, d in areal_density_df.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['areal_density'] * 1E9, ms=4, marker=mapper[g]['m'],
            color=mapper[g]['c'], alpha=0.75, linestyle='none', markeredgecolor=cor['primary_black'],
            markeredgewidth=0.5)

ax.set_ylim([0, 15])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('membrane protein density [fg / µm$^2$]', fontsize=6)
plt.savefig('./membrane_protein_density.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
for g, d in areal_density_df.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['prot_mass'] * 1E9, ms=4, marker=mapper[g]['m'],
            color=mapper[g]['c'], alpha=0.75, linestyle='none', markeredgecolor=cor['primary_black'],
            markeredgewidth=0.5)
#
ax.set_ylim([0, 120])
ax.set_xlabel('growth rate [hr$^{-1}]', fontsize=6)
ax.set_ylabel('membrane protein mass [fg / cell]', fontsize=6)
plt.savefig('./membrane_protein_mass.pdf')

# %%
phi_mem = 0.35
m_peri = 15
rho_mem = 9
alpha = 3.3

width_range = np.linspace(0.45, 1.2, 200)
pred = phi_mem * m_peri / (rho_mem * np.pi * alpha * width_range**2)

plt.plot(width_range, pred - 0.05)
for g, d in peri_prots.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    w = mass_spec_inf[(mass_spec_inf['dataset_name'] == g[0]) & (mass_spec_inf['growth_rate_hr'] == g[1]) &
                      (mass_spec_inf['condition'] == g[2]) &
                      (mass_spec_inf['interval'] == 'median') & (mass_spec_inf['quantity'] == 'mass_spec_widths')]['lower'].values[0]

    plt.plot(w, d['mass_frac'].values,
             linestyle='none', marker=mapper[g[0]]['m'],
             color=mapper[g[0]]['c'], alpha=0.5, markeredgecolor=cor['primary_black'])


plt.ylim(0, 0.15)
