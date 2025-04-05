"""
This file reads through images on a local machine and copies all of the images 
valid for the formal final analysis. Running this will not work unless you have all
raw images properly located within each folder.  
"""
#%%
import pandas as pd
import glob
import shutil 
import tqdm
import os

# Load the valid size data
data = pd.read_csv('../mass_spectrometry/valid_experimental_metadata.csv')
root = '../mass_spectrometry'

# Move over images
for g, d in tqdm.tqdm(data.groupby(['date_collected','strain', 
                          'carbon_source', 'inducer_conc', 'replicate']),
                          desc='Copying images...'):
    if g[0] == '2024-08-27':
        fnames = glob.glob(f'{root}/{g[0]}_r1/raw/images/{g[0]}_r{g[4]}_{g[2]}_{g[1]}_{g[3]}*.tif')
    else:
        fnames = glob.glob(f'{root}/{g[0]}_r1/raw/images/{g[0]}_r1_{g[1]}*{g[2]}_*_{g[3]}_rep_{g[4]}*.tif')
    for f in fnames:
        shutil.copy(f, '../../../images/') 

#%%
# Rename images from 2024-08-27 to be consistent with others.
files = glob.glob('../../../images/2024-08-27*.tif')

for f in tqdm.tqdm(files, desc='Consistent renaming...'):
    date, rep, carbon, strain, ind = f.split('/')[-1].split('_')
    rep = int(rep[1:])
    ind, suff = ind.split('-')
    ind = int(ind)
    if strain == 'relA':
        ind_type = 'dox'
    else:
        ind_type = 'IPTG'
    new_name = f'{date}_r1_{strain}_37C_{carbon}_{ind_type}_{ind}_rep_{rep}-{suff}' 
    os.rename(f, f'../../../images/{new_name}')
