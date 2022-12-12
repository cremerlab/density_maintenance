# %%
import numpy as np
import fcsparser
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import glob
cor, pal = size.viz.matplotlib_style()
dapi_files = glob.glob(
    '../../../data/flow_cytometry/2022-11-17_wildtype/samples_Sorb_DAPI*.fcs')
no_dapi_files = glob.glob(
    '../../../data/flow_cytometry/2022-11-17_wildtype/samples_Sorb_NO DAPI*.fcs')

tube_files = glob.glob(
    '../../../data/flow_cytometry/2022-11-17_wildtype/samples_Tube*.fcs')

_d, dapi = fcsparser.parse(dapi_files[0])
_n, nodapi = fcsparser.parse(no_dapi_files[0])
_t, tube = fcsparser.parse(tube_files[0])

plt.loglog(dapi['FSC-A'], dapi['V450-A'], '.', color=cor['primary_blue'])
plt.loglog(nodapi['FSC-A'], nodapi['V450-A'], '.', color=cor['primary_black'])
plt.loglog(tube['FSC-A'], tube['V450-A'], '.', color=cor['primary_red'])
# %%
fig, ax = plt.subplots(1, 1)
dx = np.sort(dapi['FSC-A'].values)
dy = np.arange(len(dx)) / len(dx)
nx = np.sort(nodapi['FSC-A'].values)
ny = np.arange(len(nx)) / len(nx)
tx = np.sort(tube['FSC-A'].values)
ty = np.arange(len(tx)) / len(nx)

ax.plot(dx, dy, '-', color=cor['primary_blue'])
ax.plot(nx, ny, '-', color=cor['primary_black'])
ax.plot(tx, ty, '-', color=cor['primary_red'])
# ax.hist(np.log10(dapi['V450-A'].values), bins=100,
#         color=cor['primary_blue'], alpha=0.5)
# ax.hist(np.log10(nodapi['V450-A'].values), bins=100,
#         color=cor['primary_black'], alpha=0.5)
# ax.hist(np.log10(tube['V450-A'].values), bins=100,
#         color=cor['primary_red'], alpha=0.5)
# ax.set_xlim([-100, 8000])
