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
mp.cpu_count()
ROOT = '../../../../data/images'
# Load images, convert to greyscale,and filter.
dirs = np.sort(glob.glob(f'{ROOT}/wildtype/2023-09-05*'))

# %%
size_df = pd.DataFrame([])
biom = []
for direc in tqdm.tqdm(dirs):
    # Get the images
    files = glob.glob(f'{direc}/*.tif')

    # Convert to grey scale
    try:
        ims = [skimage.color.rgb2gray(skimage.io.imread(f)) for f in files]
    except ValueError:
        ims = [skimage.io.imread(f) for f in files]

    # Apply tophat filtering
    print('Filtering images....')
    filt_ims = joblib.Parallel(n_jobs=mp.cpu_count()-2)(
        joblib.delayed(size.image.tophat_filter)(im) for im in tqdm.tqdm(ims))
    print('done!')

    # Iterate through each sample and segment/measure cells
    sizes = []
    splines = []
    cells = []

    # Get date info
    date, run_no = direc.split('/')[:-1].split('_r')
    for i, f in enumerate(tqdm.tqdm(files)):
        # Parse file information
        _, fname = f.split('/')[-2:]
        strain, carbon, over_expression, inducer, inducer_conc, temp, suffix = fname.split(
            '_')
        temp = float(temp[:-1])
        if date == '2022-09-29':
            ip_dist = 0.065
        else:
            ip_dist = 0.032
        print(ip_dist)
        # Process the image
        print(f'Performing segmentation of {carbon}...')
        objs, cell_df = size.image.contour_segmentation(filt_ims[i], filter=False,
                                                        area_bounds=(0.5, 15),
                                                        ecc_bounds=0.8,
                                                        solidity_bound=0.9,
                                                        perim_bounds=(
                                                            0.1, 200),
                                                        return_cells=True,
                                                        ip_dist=ip_dist,
                                                        intensity_image=ims[i])
        print('done!')
        if len(objs) == 0:
            continue
        anatomy = size.image.assign_anatomy(
            objs, cap_radius=0.75, sept_radius=0.3)
        if len(anatomy) == 0:
            continue
        biometrics = size.image.measure_biometrics(anatomy, ip_dist=ip_dist)

        if len(biometrics) == 0:
            continue
        # Assign information
        for d in [biometrics, anatomy, cell_df]:
            d['date'] = date
            d['carbon_source'] = carbon
            d['temperature_C'] = temp
            if over_expression == 'noOE':
                over_expression = 'none'
            if inducer == 'noInd':
                inducer = 'none'
            d['overexpression'] = over_expression
            d['inducer'] = inducer
            d['inducer_conc'] = inducer_conc
            d['image'] = suffix
            d['strain'] = strain
            d['run_no'] = run_no
        size_df = pd.concat([size_df, biometrics])
        biom.append(biometrics)

        sizes.append(biometrics)
        splines.append(anatomy)
        cells.append(cell_df)
    # Form dataframes for plot generation
    cell_sizes = pd.concat(sizes, sort=False)
    cell_images = pd.concat(cells, sort=False)
    cell_splines = pd.concat(splines, sort=False)

    # Save size measurements
    cell_sizes.to_csv(f'{direc}/{date}_{run_no}_sizes.csv', index=False)
    # Generate the gallerires
    print('Generating cell galleries...')
    for g, d in cell_sizes.groupby(['carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature_C']):
        idx = '_'.join([str(_g) for _g in g])
        fname = f'./output/galleries/{date}_wildtype_{idx}_gallery.png'
        suptitle = f'{date} {strain} {g[0]} OE: {g[1]} Inducer: {g[2]} {g[3]} Temp: {g[4]} C'
        cells = cell_images[(cell_images['carbon_source'] == g[0]) &
                            (cell_images['overexpression'] == g[1]) &
                            (cell_images['inducer'] == g[2]) &
                            (cell_images['inducer_conc'] == g[3]) &
                            (cell_images['temperature_C'] == g[4])]
        splines = cell_splines[(cell_splines['carbon_source'] == g[0]) &
                               (cell_splines['overexpression'] == g[1]) &
                               (cell_splines['inducer'] == g[2]) &
                               (cell_splines['inducer_conc'] == g[3]) &
                               (cell_splines['temperature_C'] == g[4])]
        _ = size.viz.cell_gallery(
            d, cells, splines, fname=fname, suptitle=suptitle)
        plt.close()
    print('done!')


# %%

# %%
# Save the huge size dataframe
size_df = pd.concat([pd.read_csv(f)
                    for f in glob.glob(f'{ROOT}/*/*/*.csv')], sort=False)
size_df.to_csv('./output/compiled_size_measurements.csv')
