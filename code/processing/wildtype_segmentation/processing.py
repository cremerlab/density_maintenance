#%%
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
files = np.sort(glob.glob('../../../data/images/wt_test/*.tif'))
ims = [skimage.color.rgb2gray(skimage.io.imread(f)) for f in files]
#%%
mp.cpu_count()
# %%
carbon_dict = {'Acetate':'acetate',
               'glu': 'glucose',
               'Glu_CAA':'glucose + CAA',
               'glycrk': 'glycerol',
               'LB': 'LB'}

filt_ims = joblib.Parallel(n_jobs=mp.cpu_count()-1)(
                    joblib.delayed(size.image.tophat_filter)(im) for im in tqdm.tqdm(ims))

#%%
sizes = []
splines = []
for i, f in enumerate(tqdm.tqdm(files)):
    # Parse file information
    carbon = f.split('/')[-1].split('WT')[0]
    image = f.split('/')[-1].split('-')[-1].split('.tif')[0]

    # Process the image 
    objs, cells = size.image.contour_segmentation(filt_ims[i], filter=False, 
                                    area_bounds=(0.5, 10), ecc_bounds=0.5,
                                    solidity_bound=0.9,  return_cells=True,
                                    intensity_image=ims[i])
    if len(objs) == 0:
        continue
    anatomy = size.image.assign_anatomy(objs, max_curve=1)
    biometrics = size.image.measure_biometrics(anatomy)

    # Assign information
    biometrics['carbon_source'] = carbon_dict[carbon]
    biometrics['image'] = image
    biometrics['strain'] = 'WT'
    biometrics['image'] = image
    biometrics['cell_image'] = [cells[idx]['intensity_image'] for idx in biometrics['cell_id'].values]
    anatomy['carbon_source'] = carbon_dict[carbon]
    anatomy['image'] = image
    anatomy['strain'] = 'WT'
    sizes.append(biometrics)
    splines.append(anatomy)

cell_sizes = pd.concat(sizes, sort=False)

#%%    
cell_sizes['idx'] = biometrics.groupby(['carbon_source', 'cell_id', 'image']).ngroup() + 1

n_cells = biometrics['idx'].max()
cols = 5
rows = int(np.ceil(n_cells/cols))

fig, ax = plt.subplots(rows, cols, sharex=True, sharey=True)
ax = ax.ravel()
for a in ax:
    a.axis('off')
    a.set_aspect('equal')

for g, d in biometrics.groupby(['cell_id', 'image', 'idx']):
    cell = anatomy[(anatomy['cell_id']==g[0]) & (anatomy['image']==g[1])]
    ax[g[2]-1].imshow(d['cell_image'].values[0], cmap='Greys_r')
    sides = cell[(cell['component']=='left') | (cell['component']=='right')] 
    caps = cell[(cell['component']=='top') | (cell['component']=='bottom')] 
    ax[g[2]-1].plot(sides['x_coords'], sides['y_coords'], color=cor['primary_green'], marker='.',
                markeredgewidth=0, lw=0, ms=2)
    ax[g[2]-1].plot(caps['x_coords'], caps['y_coords'], color=cor['primary_blue'], marker='.',
                markeredgewidth=0, lw=0, ms=2)

#%% 
# %%
for g, d in cell_sizes.groupby(['carbon_source']):
    x = np.sort(d['width_mean'])
    y = np.arange(len(x))/len(x)
    plt.plot(x, y, label=g)
plt.xlabel('average cell width [µm]')
plt.ylabel('cumulative distribution')
plt.legend()
plt.savefig('./cell_width_distributions.pdf')
# %%
