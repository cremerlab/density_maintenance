# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')

fig, ax = plt.subplots(1, 3, figsize=(4.5, 1.5), sharex=True)
for a in ax:
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
    a.set_xlim([0, 2.1])
ax[0].set_ylabel('w [µm]\naverage cell width', fontsize=6)
ax[0].set_ylim([0.45, 1.2])
ax[1].set_ylim([1, 3.5])
ax[2].set_ylim([2, 9])
ax[1].set_ylabel('$\ell$ [µm]\naverage cell length', fontsize=6)
ax[2].set_ylabel('$S_A/V$ [µm]\nsurface-to-volume', fontsize=6)
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['width_um'], **fmt, ms=4)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], **fmt, ms=4)
    ax[2].plot(d['growth_rate_hr'], d['surface_to_volume'], **fmt, ms=4)

plt.tight_layout()
ax[0].legend(bbox_to_anchor=(3.7, -0.35), fontsize=6, ncols=3, frameon=False)
fig.text(0.01, 0.95, '(A)', fontsize=8)
fig.text(0.33, 0.95, '(B)', fontsize=8)
fig.text(0.66, 0.95, '(C)', fontsize=8)
plt.savefig('../../figures/FigS1_size_relations.pdf', bbox_inches='tight')
