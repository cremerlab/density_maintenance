# %%
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.analytical
import seaborn as sns
cor, pal = size.viz.matplotlib_style()
cmap = sns.color_palette("ch:start=.2,rot=-.3", n_colors=8)

data = pd.read_csv('../../data/mcmc/wildtype_hyperparameter_size_samples.csv')
# data['width_um'] -= 0.3
voldata = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
grdata = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')
# %%
_dfs = []
data['periplasmic_volume_fL'] = size.analytical.envelope_volume(
    data['length_um'], data['width_um'], 0.025)

# Define the percentiles to compute for each hyperparameter
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']

growth_rates = {'acetate': {'growth_rate': 0.35,
                            'protein': 23},
                'sorbitol': {'growth_rate': 0.46,
                             'protein': 18},
                'glycerol': {'growth_rate': 0.68,
                             'protein': 16.8},
                'glucose': {'growth_rate': 0.87,
                            'protein': 10.76},
                'glucoseCAA': {'growth_rate': 1.17,
                               'protein': 10.33}}


def cell_number(vol):
    return 1/(1.3372E-9 * vol - 4.71E-10)


for g, d in data.groupby(['carbon_source']):
    if g in growth_rates.keys():
        # Compute various summary statistics of the hyper parameters
        prot_per_cell = growth_rates[g]['protein'] / \
            cell_number(d['volume_fL']) * 1E9
        density = prot_per_cell / d['periplasmic_volume_fL']
        _percs = np.percentile(density, percs)
        _perc_df = pd.DataFrame([_percs], columns=perc_cols)
        _perc_df['mean'] = density.mean()
        _perc_df['median'] = density.median()
        _perc_df['carbon_source'] = g
        _perc_df['strain'] = 'wildtype'
        _dfs.append(_perc_df)
df = pd.concat(_dfs, sort=False)

# %%
# Do a quick eexponential regression on the size medians to properly normalize mass spec data
vol, gr, ell, w = [], [], [], []
for g, d in voldata.groupby(['carbon_source']):
    _vol = voldata[(voldata['carbon_source'] == g) &
                   (voldata['parameter'] == 'volume_fL')]
    _ell = voldata[(voldata['carbon_source'] == g) &
                   (voldata['parameter'] == 'length_um')]
    _w = voldata[(voldata['carbon_source'] == g) &
                 (voldata['parameter'] == 'width_um')]
    _gr = grdata[(grdata['carbon_source'] == g)]

    if len(_gr) != 0:
        vol.append(_vol['median'].values[0])
        gr.append(_gr['median'].values[0])
        ell.append(_ell['median'].values[0])
        w.append(_w['median'].values[0])
vol_popt = scipy.stats.linregress(gr, np.log(vol))
ell_popt = scipy.stats.linregress(gr, np.log(ell))
w_popt = scipy.stats.linregress(gr, np.log(w))


def lam2val(lam, popt):
    return np.exp(popt[0]*lam + popt[1])


# %%
mass_spec = pd.read_csv(
    '../../data/calculated_envelope_protein_mass_and_concentration.csv')
mass_spec['length'] = lam2val(mass_spec['growth_rate_hr'], ell_popt)
mass_spec['width'] = lam2val(mass_spec['growth_rate_hr'], w_popt)
mass_spec['volume'] = lam2val(mass_spec['growth_rate_hr'], vol_popt)
mass_spec['V_peri'] = size.analytical.envelope_volume(
    mass_spec['length'], mass_spec['width'], delta=0.025)
mass_spec['phi_peri'] = mass_spec['periplasmic_protein_mass_fg'].values / \
    mass_spec['total_protein_mass_fg'].values
mass_spec['tot_prot'] = 0.55 * 0.3 * 1.1 * mass_spec['volume']
mass_spec['m_peri'] = mass_spec['phi_peri'] * mass_spec['tot_prot']
mass_spec['rho_peri'] = (mass_spec['m_peri'] / mass_spec['V_peri']) * 1E3

# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 4))
gr_range = np.linspace(0.3, 1.3, 200)
ax[0].plot(gr_range, lam2val(gr_range, vol_popt), '--')
ax[1].plot(mass_spec['growth_rate_hr'], mass_spec['rho_peri'], 'o')
mapper = {'acetate': cmap[6],
          'sorbitol': cmap[5],
          'glycerol': cmap[4],
          'glucose': cmap[3],
          'glucoseCAA': cmap[2],
          'LB': cmap[0],
          'ezMOPS': cmap[1]}
for g, d in df.groupby(['carbon_source']):
    _vol = voldata[(voldata['carbon_source'] == g) &
                   (voldata['parameter'] == 'volume_fL')]
    _gr = grdata[grdata['carbon_source'] == g]
    ax[1].vlines(_gr['median'], d['2.5%'],
                 d['97.5%'], lw=0.75, color=mapper[g])
    ax[1].vlines(_gr['median'], d['12.5%'],
                 d['87.5%'], lw=1.5, color=mapper[g])
    ax[1].hlines(d['median'], _gr['2.5%'],
                 _gr['97.5%'], lw=0.75, color=mapper[g])
    ax[1].hlines(d['median'], _gr['12.5%'],
                 _gr['87.5%'], lw=1.5, color=mapper[g])
    ax[0].vlines(_gr['median'], _vol['2.5%'], _vol['97.5%'],
                 lw=0.75, color=mapper[g])
    ax[0].vlines(_gr['median'], _vol['12.5%'], _vol['87.5%'],
                 lw=1.5, color=mapper[g])
    ax[0].hlines(_vol['median'], _gr['2.5%'], _gr['97.5%'],
                 lw=0.75, color=mapper[g])
    ax[0].hlines(_vol['median'], _gr['12.5%'], _gr['87.5%'],
                 lw=1.5, color=mapper[g])
    ax[0].plot(_gr['median'], _vol['median'], marker='o', ms=3, markerfacecolor='white',
               markeredgecolor=mapper[g], markeredgewidth=1.5)
    ax[1].plot(_gr['median'], d['median'], marker='o', ms=3, markerfacecolor='white',
               markeredgecolor=mapper[g], markeredgewidth=1.5)

ax[1].set_ylim([0, 200])
# %%
fig, ax = plt.subplots(1, 1)
for g, d in data.groupby(['carbon_source']):
    ax.hist(d['width_um'], bins=25, alpha=0.5, label=g)
ax.legend()

# %%
for g, d in data.groupby(['carbon_source']):
    plt.hist(d['width_um'], bins=20, alpha=0.5, label=g)
plt.legend()

# %%
size_data = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')
glyc = size_data[(size_data['carbon_source'] == 'glycerol')]  # &
#  (size_data['date'] == '2022-03-11_')]
for g, d in glyc.groupby(['date']):
    plt.hist(d['width_median'], alpha=0.5, bins=20, label=g)
plt.legend()
