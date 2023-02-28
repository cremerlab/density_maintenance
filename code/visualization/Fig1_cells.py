# %%
import numpy as np
import pandas as pd
import skimage.io
import skimage.color
import size.image
import matplotlib.pyplot as plt
import tqdm
import size.viz
import glob
cor, pal = size.viz.matplotlib_style()


# Define the data directories, using only 1
ac_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*acetate*.tif')[0]
glu_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*glucose*.tif')[0]
sorb_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*sorbitol*.tif')[0]
gly_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*glycerol*.tif')[0]
lb_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*LB*.tif')[0]
files = np.array([ac_files, glu_files, sorb_files,
                 gly_files, lb_files]).flatten()

# Iterate
cell_df = pd.DataFrame([])
anat_df = pd.DataFrame([])
biom_df = pd.DataFrame([])
for i, f in enumerate(tqdm.tqdm(files)):
    carb = f.split('/')[-1].split('_')[1]
    im = skimage.io.imread(f)
    im = skimage.color.rgb2gray(im)
    obj, cells = size.image.contour_segmentation(im, return_cells=True)
    anatomy = size.image.assign_anatomy(obj, cap_radius=1)
    biometrics = size.image.measure_biometrics(anatomy)
    cells['carbon'] = carb
    cells['image'] = i
    anatomy['carbon'] = carb
    anatomy['image'] = i
    biometrics['carbon'] = carb
    biometrics['image'] = i
    biom_df = pd.concat([biom_df, biometrics])
    cell_df = pd.concat([cell_df, cells])
    anat_df = pd.concat([anat_df, anatomy])

# %%
fig, ax = plt.subplots(4, 3, figsize=(4, 2), sharex=True, sharey=True)

for a in ax.ravel():
    a.set_xticks([])
    a.set_yticks([])
    a.grid(False)
    a.set_xlim(0, 200)
    a.set_facecolor('none')

ax[0, 0].set_title('acetate', fontsize=6)
ax[0, 1].set_title('glucose', fontsize=6)
ax[0, 2].set_title('LB', fontsize=6)

col = {'acetate': 0, 'glucose': 1, 'LB': 2}
for g, d in biom_df.groupby(['carbon']):
    if g not in ['acetate', 'glucose', 'LB']:
        continue
    d = d.sort_values(by='length')
    inds = [0, int(len(d)/4),  int(3 * len(d)/4), -1]
    for i, j in enumerate(inds):

        # Isolate the cell
        cell = cell_df[(cell_df['carbon'] == g) &
                       (cell_df['cell_id'] == d['cell_id'].values[j])]['cell_image'].values[0]
        cell = (cell - cell.min()) / (cell.max() - cell.min())
        # Isolate the anatomy
        anat = anat_df[(anat_df['carbon'] == g) &
                       (anat_df['cell_id'] == d['cell_id'].values[j])]

        # Plot the cell and the contour
        ax[i, col[g]].imshow(cell.T, cmap='Greys_r', vmin=0, vmax=1.1)

        for _g, _d in anat.groupby(['component']):
            if _g in ['top', 'bottom']:
                c = cor['purple']
            else:
                c = cor['pale_blue']
            ax[i, col[g]].plot(_d['y_coords'][::50], _d['x_coords'][::50], '.', markeredgewidth=0, color=c,
                               ms=1.5)


plt.subplots_adjust(hspace=-0.7)
plt.savefig('../../figures/Fig1_contours.pdf')
