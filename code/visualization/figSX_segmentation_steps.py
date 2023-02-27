# %%

import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import scipy.ndimage
import scipy.interpolate
import skimage.morphology
import skimage.color
import size.image
import size.viz
cor, pal = size.viz.matplotlib_style()

# Load image and crop to region to demonstrate segmentation
image = skimage.io.imread(
    '../../data/example_images/glucose2.tif')
image = skimage.color.rgb2gray(image[:, :, :-1])
plt.imshow(image)
# %%
# Save naked image for display
plt.imshow(image, cmap='Greys_r')
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_01_phase.pdf',
            bbox_inches='tight')

# %%
# Filtering step
im_filt = size.image.tophat_filter(image, threshold='none')
plt.imshow(im_filt)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_02_tophat.pdf')

# %%
# Thresholding step
im_filt = size.image.tophat_filter(image)
plt.imshow(im_filt)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_03_otsu.pdf')

# %%
# LOG segmentation
selem = skimage.morphology.square(2)
seg = size.image.log_segmentation(im_filt, selem=selem,
                                  radius=2, thresh=0.0001, median_filt=False,
                                  label=True)
seg = skimage.morphology.remove_small_holes(seg)
seg = skimage.morphology.remove_small_objects(seg)
labeled = skimage.measure.label(seg)
# labeled[labeled == 2] = 0
plt.imshow(labeled)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_04_LoG.pdf')

# %%
# Cell isolation
plt.imshow(labeled == 6, cmap='Greys_r')
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_05_isolated.pdf')

# %%
# Cell rotation
props = skimage.measure.regionprops(labeled)
cell = props[4]
padded, _ = size.image.pad_bbox(cell.bbox, np.shape(labeled == 6), pad=10)
rot = scipy.ndimage.rotate(
    seg[padded], -np.rad2deg(props[4].orientation), order=0) > 0
plt.imshow(rot, cmap='Greys_r')
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_06_rotated.pdf')


# %%
# Morphological operations
erode_selem = skimage.morphology.rectangle(1, 5)
rot = skimage.morphology.binary_erosion(rot, erode_selem)
rot = skimage.morphology.remove_small_holes(rot)
rot = skimage.morphology.remove_small_objects(rot)
rot = scipy.ndimage.binary_fill_holes(rot)
relab = skimage.measure.label(rot.astype(int))
plt.imshow(relab, cmap='Greys_r')
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_07_morph.pdf')

# %%
# Contouring
cell_grey = scipy.ndimage.rotate(
    image[padded], -np.rad2deg(props[4].orientation), order=5, mode='nearest')
rot_props = skimage.measure.regionprops(relab)
bbox = rot_props[0].bbox
rot_pad, _ = size.image.pad_bbox(bbox, np.shape(rot), pad=10)
cont = skimage.measure.find_contours(rot[rot_pad], 0)[0]
plt.imshow(rot[rot_pad], cmap='Greys_r')
plt.plot(cont[:, 1], cont[:, 0], '.', ms=8, color=cor['primary_red'], lw=3)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_08_cont.pdf')

# %%
# Contour filtering
cx = scipy.ndimage.uniform_filter(cont[:, 1], 10, mode='wrap')
cy = scipy.ndimage.uniform_filter(cont[:, 0], 10, mode='wrap')
plt.imshow(cell_grey[rot_pad], cmap='Greys_r')
plt.plot(cx, cy, '.', ms=8, color=cor['primary_red'], lw=3)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_09_contfilt.pdf')

# %%
# Interpolation
tck, _ = scipy.interpolate.splprep([cx, cy], per=1, k=3, s=30)
unew = np.arange(0, 1.0001, 0.001)
out = scipy.interpolate.splev(unew, tck)
x = out[0][:-2]
y = out[1][:-2]
plt.imshow(cell_grey[rot_pad], cmap='Greys_r')
plt.plot(x, y, '-', ms=8, color=cor['primary_red'], lw=3)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_10_spline.pdf')
# %%
# Compute curvature
k = size.image.compute_curvature(out)
plt.imshow(cell_grey[rot_pad], cmap='Greys_r')
plt.scatter(x, y, c=k/0.065, s=30)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_11_curvature.pdf')

# %%
out = size.image.contour_segmentation(image, ip_dist=0.035)

# %%
# Assign domains based on curvature
cell_cont = out[out['cell_id'] == 4]
anat = size.image.assign_anatomy(cell_cont, cap_radius=1)
plt.imshow(cell_grey[rot_pad], cmap='Greys_r')
for g, d in anat.groupby(['component']):
    if g in ['top', 'bottom']:
        c = cor['primary_blue']
    else:
        c = cor['primary_green']
    plt.plot(d['x_coords'], d['y_coords'], '.',
             ms=5, markeredgewidth=0, color=c)
plt.xticks([])
plt.yticks([])
plt.savefig('../../figures/supplement/seg_dem_12_anatomy.pdf')
