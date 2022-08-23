# %%
import numpy as np
import pandas as pd
import skimage.color
import matplotlib.pyplot as plt
import size.image
import glob
import tqdm
import size.viz
import joblib
import multiprocessing as mp
cor, pal = size.viz.matplotlib_style()

DATE = ""
STRAIN = 'WT'

files = np.sort(glob.glob('../../../data/images/wt_test/*.tif'))
ims = [skimage.color.rgb2gray(skimage.io.imread(f)) for f in files]
# %%
mp.cpu_count()
carbon_dict = {'Acetate': 'acetate',
               'glu': 'glucose',
               'Glu_CAA': 'glucose + CAA',
               'glycrk': 'glycerol',
               'LB': 'LB'}
filt_ims = joblib.Parallel(n_jobs=mp.cpu_count()-1)(
    joblib.delayed(size.image.tophat_filter)(im) for im in tqdm.tqdm(ims))

# %%
sizes = []
splines = []
cells = []
for i, f in enumerate(tqdm.tqdm(files)):
    # Parse file information
    carbon = f.split('/')[-1].split('WT')[0]
    image = f.split('/')[-1].split('-')[-1].split('.tif')[0]

    # Process the image
    objs, cell_df = size.image.contour_segmentation(filt_ims[i], filter=False,
                                                    area_bounds=(0.5, 10), ecc_bounds=0.8,
                                                    solidity_bound=0.9,  return_cells=True,
                                                    intensity_image=ims[i])
    if len(objs) == 0:
        continue
    anatomy = size.image.assign_anatomy(objs, cap_radius=1)
    biometrics = size.image.measure_biometrics(anatomy)

    if len(biometrics) == 0:
        continue
    # Assign information
    for d in [biometrics, anatomy, cell_df]:
        d['carbon_source'] = carbon_dict[carbon]
        d['image'] = image
        d['strain'] = 'WT'

    sizes.append(biometrics)
    splines.append(anatomy)
    cells.append(cell_df)
cell_sizes = pd.concat(sizes, sort=False)
cell_images = pd.concat(cells, sort=False)
cell_splines = pd.concat(splines, sort=False)
cell_sizes.to_csv('../../../wt_cell_sizes.csv', index=False)
# %%

for g, d in cell_sizes.groupby(['carbon_source']):
    fname = f'./output/{g}_gallery.png'
    suptitle = f'segmented cells – growth on {g}'
    cells = cell_images[cell_images['carbon_source'] == g]
    splines = cell_splines[cell_splines['carbon_source'] == g]
    _ = size.viz.cell_gallery(
        d, cells, splines, fname=fname, suptitle=suptitle)

# %%
