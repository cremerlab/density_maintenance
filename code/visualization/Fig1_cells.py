# %%
import numpy as np
import pandas as pd
import skimage.io
import skimage.color
import size.image
import tqdm
import size.viz
import glob
cor, pal = size.viz.matplotlib_style()


# Define the data directories, using only 5 images
ac_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*acetate*.tif')[:5]
glu_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*glucose*.tif')[:5]
sorb_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*sorbitol*.tif')[:5]
gly_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*glycerol*.tif')[:5]
lb_files = glob.glob(
    '../../data/images/wildtype/2022-04-14_r1/*LB*.tif')[:5]
files = np.array([ac_files, glu_files, sorb_files,
                 gly_files, lb_files]).flatten()

# Iterate
cell_df = pd.DataFrame([])
anat_df = pd.DataFrame([])
for i, f in enumerate(tqdm.tqdm(files)):
    carb = f.split('/')[-1].split('_')[1]
    for j, s in enumerate(tqdm.tqdm(f)):
        im = skimage.io.imread(f)
        im = skimage.color.rgb2gray(im)
        obj, cells = size.image.contour_segmentation(im, return_cells=True)
        anatomy = size.image.assign_anatomy(obj, cap_radius=1)
        cells['carbon'] = carb
        cells['image'] = j
        anatomy['carbon'] = carb 
        anatomy['image'] = j
        cell_df = pd.concat([cell_df, cells])
        anat_df = pd.concat([anat_df, anatomy])

#%%

