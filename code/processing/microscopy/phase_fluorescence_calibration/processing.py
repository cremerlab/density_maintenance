# %%
import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.image
import size.viz
import skimage.io
import scipy.ndimage
import glob

import tqdm
cor, pal = size.viz.matplotlib_style()

# Define the phase and fluorescence files
phase_ims = skimage.io.ImageCollection(list(np.sort(glob.glob(
    '../../../../data/images/phase_fluorescence_calibration/*phase*.tif'))))
fluo_ims = skimage.io.ImageCollection(list(np.sort(glob.glob(
    '../../../../data/images/phase_fluorescence_calibration/*fluorescence*.tif'))))

# Filter the phase images
plt.imshow(fluo_ims[-3])
plt.imshow(phase_ims[-3])

# %%
imp.reload(size.image)
coords, mask, cells = size.image.contour_segmentation(
    phase_ims[-3], return_mask=True, return_cells=True, ip_dist=0.065,
    area_bounds=(1, 100))  # , intensity_image=fluo_ims[-3])

# %%
anatomy = size.image.assign_anatomy(coords, cap_radius=1)
biometrics = size.image.measure_biometrics(anatomy, ip_dist=0.065)
biometrics['image'] = 1
anatomy['image'] = 1
cells['image'] = 1


# %%
gal = size.viz.cell_gallery(biometrics, cells, anatomy, './test.png')

# %%
coords, mask, cells = size.image.contour_segmentation(
    fluo_ims[-3], filter=False, return_mask=True, return_cells=True, ip_dist=0.065,
    area_bounds=(0.1, 1000), ecc_bound=0.9, solidity_bound=0.5)
plt.imshow(mask)

# %%
anatomy = size.image.assign_anatomy(coords, cap_radius=5)
biometrics = size.image.measure_biometrics(anatomy, ip_dist=0.065)
biometrics['image'] = 1
anatomy['image'] = 1
cells['image'] = 1

# %%
gal = size.viz.cell_gallery(biometrics, cells, anatomy, './test_fluo.pdf')

# %%
im_float = skimage.img_as_float(fluo_ims[-3])
im_LoG = scipy.ndimage.filters.gaussian_laplace(im_float, 2)

# Try a different way of
selem = skimage.morphology.square(3)
im_max = scipy.ndimage.filters.maximum_filter(im_LoG, footprint=selem)
im_min = scipy.ndimage.filters.minimum_filter(im_LoG, footprint=selem)
sobel = skimage.filters.sobel(im_LoG)
zero_cross = (((im_LoG >= 0) & (im_min < 0)) | ((im_LoG <= 0) & (im_max > 0)))\
    & (sobel >= 0.0001)
skel = skimage.morphology.skeletonize(zero_cross)
plt.imshow(skel)


# %%

# %%
