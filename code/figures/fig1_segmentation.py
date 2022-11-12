# %%
import numpy as np
import pandas as pd
import skimage.io
import size.image
import size.viz
import glob
import matplotlib.pyplot as plt
import tqdm
cor, pal = size.viz.matplotlib_style()
carbs = ['acetate', 'glucose', 'LB']
anat_dfs = []
phase_cell_dfs, fluo_cell_dfs = [], []
for i, c in enumerate(tqdm.tqdm(carbs)):
    # Load some example images.
    phase = skimage.io.ImageCollection(
        glob.glob(f'../../data/images/phase_fluorescence_calibration/*{c}*phase*'))
    fluo = skimage.io.ImageCollection(
        glob.glob(f'../../data/images/phase_fluorescence_calibration/*{c}*fluo*'))
    coords, mask, phase_cells = size.image.contour_segmentation(
        phase[-1], return_mask=True, return_cells=True, ip_dist=0.065,
        area_bounds=(1, 100))
    _, _, fluo_cells = size.image.contour_segmentation(
        phase[-1], return_mask=True, return_cells=True, ip_dist=0.065,
        area_bounds=(1, 100), intensity_image=fluo[-1])

    anatomy = size.image.assign_anatomy(coords, cap_radius=1)
    biometrics = size.image.measure_biometrics(anatomy, ip_dist=0.065)
    biometrics['image'] = 1
    anatomy['image'] = 1
    phase_cells['image'] = 1
    fluo_cells['image'] = 1
    anat_dfs.append(anatomy)
    phase_cell_dfs.append(phase_cells)
    fluo_cell_dfs.append(fluo_cells)

# %%
fig, ax = plt.subplots(2, 3, figsize=(3, 3))  # , sharex=True, sharey=True)
for a in ax.ravel():
    a.axis('off')
cell_ids = [1, 15, 1]
for i in range(3):
    _phase = phase_cell_dfs[i]
    _fluo = fluo_cell_dfs[i]
    _anat = anat_dfs[i][anat_dfs[i]['cell_id'] == cell_ids[i]]
    ax[0, i].imshow(_phase[_phase['cell_id'] == cell_ids[i]]['cell_image'].values[0],
                    cmap='Greys_r')
    for g, d in _anat.groupby(['component']):
        ax[0, i].plot(d['x_coords'].values[::10], d['y_coords'].values[::10], '.',
                      markeredgewidth=0, ms=1, color=cor['primary_red'], lw=1)
        ax[1, i].plot(d['x_coords'].values[::10], d['y_coords'].values[::10], '.',
                      markeredgewidth=0, ms=1, color=cor['primary_red'], lw=1)

    ax[1, i].imshow(_fluo[_fluo['cell_id'] == cell_ids[i]]['cell_image'].values[0],
                    cmap='mako')

for a in ax.ravel():
    a.set_xlim([0, 3 / 0.065])
# phase_cell_dfs
plt.savefig('./fig1_segmentation.pdf')
