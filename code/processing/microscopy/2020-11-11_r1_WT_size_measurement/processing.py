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

ROOT = '../../../../data/images/'
RUN_NO = 1
DATE = "2020-11-11"
STRAIN = 'WT'
carbon_dict = {'Acetate': 'acetate',
               'glu': 'glucose',
               'Glu_CAA': 'glucose + CAA',
               'glycrk': 'glycerol',
               'LB': 'LB'}

# %%
# Load images, convert to greyscale,and filter.
files = np.sort(
    glob.glob(f'../../../../data/images/{DATE}_r{RUN_NO}_{STRAIN}/*.tif'))
ims = [skimage.color.rgb2gray(skimage.io.imread(f)) for f in files]
mp.cpu_count()
filt_ims = joblib.Parallel(n_jobs=mp.cpu_count()-2)(
    joblib.delayed(size.image.tophat_filter)(im) for im in tqdm.tqdm(ims))

# %%
# Iterate through each sample and segment/measure cells
sizes = []
splines = []
cells = []
for i, f in enumerate(tqdm.tqdm(files)):
    # Parse file information
    carbon = f.split('/')[-1].split(f'{STRAIN}')[0]
    image = f.split('/')[-1].split('-')[-1].split('.tif')[0]

    # Process the image
    objs, cell_df = size.image.contour_segmentation(filt_ims[i], filter=False,
                                                    area_bounds=(0.5, 10),
                                                    ecc_bounds=0.8,
                                                    solidity_bound=0.9,
                                                    return_cells=True,
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
        d['strain'] = STRAIN
        d['date'] = DATE

    sizes.append(biometrics)
    splines.append(anatomy)
    cells.append(cell_df)
cell_sizes = pd.concat(sizes, sort=False)
cell_images = pd.concat(cells, sort=False)
cell_splines = pd.concat(splines, sort=False)

# %%
cell_sizes.to_csv(
    f'./output/{DATE}_r{RUN_NO}_{STRAIN}_cell_sizes.csv', index=False)
cell_splines.to_csv(
    f'./output/{DATE}_r{RUN_NO}_{STRAIN}_cell_splines.csv', index=False)
# %%
# Generate summary plots
for g, d in cell_sizes.groupby(['carbon_source']):
    fname = f'./output/{DATE}_r{RUN_NO}_{STRAIN}_{g}_gallery.png'
    suptitle = f'{DATE} – run {RUN_NO} – {STRAIN} growth on {g}'
    cells = cell_images[cell_images['carbon_source'] == g]
    splines = cell_splines[cell_splines['carbon_source'] == g]
    _ = size.viz.cell_gallery(
        d, cells, splines, fname=fname, suptitle=suptitle)

# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 2))
for a in ax:
    a.set_ylabel('cumulative distribution')
ax[0].set_xlabel('cell width [µm]')
ax[1].set_xlabel('cell length [µm]')

for g, d in cell_sizes.groupby('carbon_source'):
    y = np.arange(len(d))/len(d)
    width = np.sort(d['width_median'].values)
    length = np.sort(d['length'].values)
    ax[0].plot(width, y, '-', lw=1, label=g)
    ax[1].plot(length, y, '-', lw=1, label=g)
ax[0].legend()
plt.savefig(f'./output/{DATE}_r{RUN_NO}_{STRAIN}_size_distributions.png')

# %%
