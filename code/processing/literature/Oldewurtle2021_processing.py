# %%
import numpy as np
import pandas as pd
import scipy.io
import glob
files = glob.glob('../../../data/literature/Oldewurtle2021/*.mat')

mat = scipy.io.loadmat(files[0], squeeze_me=True)
np.shape(mat['roiData'])
