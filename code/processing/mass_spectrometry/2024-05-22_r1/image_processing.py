#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import joblib 
import multiprocessing as mp
import tqdm
import skimage.color
import size.viz 
import size.image
import glob
cor, pal = size.viz.matplotlib_style()

IP_DIST = 0.032
DATE = "2024-05-22"
CARBON_SOURCES = ['glucose', 'glucose+acetate',
                  'glycerol', 'sorbitol']
size_df = pd.DataFrame([])
for i in range(1, 3):
    for j, carb in enumerate(CARBON_SOURCES):
        files = glob.glob(f'./raw/images/*{carb}*_rep_{i}-*.tif')
        if len(files) == 0:
            continue
        # Load and cast the images to grayscale
        ims = [skimage.color.rgb2gray(skimage.io.imread(f)) for f in files]

        # Filter the images  
        filt_ims = joblib.Parallel(n_jobs=mp.cpu_count()-2)(
            joblib.delayed(size.image.tophat_filter)(im) for im in tqdm.tqdm(ims)) 
    
        # Iterate through each and fit splines and cells
        sizes, cells, splines = [], [], []
        for k, f in enumerate(tqdm.tqdm(files)):
            # Unpack the image name
            date, run_no, strain, temp_C, carbon, ind, ind_conc, _, rep_suff = f.split('/')[-1].split('_')
            temp_C = float(temp_C[:-1])
            run_no = float(run_no[1])
            rep, suff = rep_suff.split('-')
            suff = suff.split('.')[0]

            # Segment
            objs, cell_df = size.image.contour_segmentation(filt_ims[k], filter=False,
                                                        area_bounds=(0.5, 15),
                                                        ecc_bound=0.8,
                                                        solidity_bound=0.9,
                                                        perim_bounds=(0, 100),
                                                        return_cells=True,
                                                        ip_dist=IP_DIST,
                                                        intensity_image=ims[k])
            # If segmentation worked, try to assign anatomy
            if len(objs) == 0:
                continue
            anatomy = size.image.assign_anatomy(objs, cap_radius=0.75, sept_radius=0.3)

            # If anatomy assignment worked, try to measure biometrics
            if len(anatomy) == 0:
                continue 
            biometrics = size.image.measure_biometrics(anatomy, ip_dist=IP_DIST)

            # If biometrics worked, save the data
            if len(biometrics) == 0:
                continue

            for d in [biometrics, anatomy, cell_df]:
                d['date'] = date
                d['replicate'] = rep
                d['run_no'] = run_no
                d['strain'] = strain
                d['temperature_C'] = temp_C
                d['carbon_source'] = carbon
                d['inducer'] = ind
                d['inducer_conc'] = ind_conc
                d['image'] = suff
                sizes.append(biometrics)
                splines.append(anatomy)
                cells.append(cell_df)

        # Form the dataframes for the plot generation
        cell_sizes = pd.concat(sizes, sort=False)
        cell_images = pd.concat(cells, sort=False)
        cell_splines = pd.concat(splines, sort=False)

        # Add the cell size info to the size dataframe
        size_df = pd.concat([size_df, cell_sizes], sort=False)

        # Generate the galleries
        fname = f'./viz/{DATE}_{carb}_rep{i}_gallery.png'
        suptitle=f'{strain} {carbon} replicate {rep}'
        _ = size.viz.cell_gallery(cell_sizes, cell_images, cell_splines,
                              fname=fname, suptitle=suptitle)
        plt.close()
        
size_df.to_csv(f'./processed/{DATE}_r1_sizes.csv', index=False)

#%%
# Generate CDFs
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.set_ylabel('empirical CDF', fontsize=6)
for g, d in size_df.groupby(['carbon_source','replicate']):
    y = np.arange(len(d))/len(d)
    ax[0, 0].plot(np.sort(d['width_median'].values), y, '-', lw=1, label=g)
    ax[0, 1].plot(np.sort(d['length'].values), y, '-', lw=1, label=g)
    ax[1, 0].plot(np.sort(d['surface_to_volume'].values), y, '-', lw=1, label=g)
    ax[1, 1].plot(np.sort(d['volume'].values), y, '-', lw=1, label=g)
ax[0,0].set_xlabel('width [µm]', fontsize=6)
ax[0,1].set_xlabel('length [µm]', fontsize=6)
ax[1,0].set_xlabel('SA/V [µm$^{-1}$]', fontsize=6)
ax[1,1].set_xlabel('volume [µm$^3$]', fontsize=6)
leg = ax[0, 0].legend(title='replicate')
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig(f'./viz/{DATE}_r1_size_cdfs.png')