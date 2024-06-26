#
# Loop through all the DISORT output and reduce its dimensionality by picking a slope and aspect
#

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interpn
from pyRT_DISORT_MCD_utils import MakeFluxLookupTable

#
# Setup Ls values and all the filenames to consolidate into four tables
#
lsubs    = np.arange(0, 361, 15)
fnames   = ['../pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_' + f"{item:03}" + '.pkl' for item in lsubs]
vnames   = ['../pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_' + f"{item:03}" + '.pkl' for item in lsubs]

slope    = np.deg2rad(70.0)
aspect   = np.deg2rad(225.0)
for i,savename in enumerate(fnames):        
    direct, diffuse, atmless, direct_axes, diffuse_axes = MakeFluxLookupTable(savename, slope, aspect)

slope    = np.deg2rad(00.0)
aspect   = np.deg2rad(225.0)
for i,savename in enumerate(fnames):        
    direct, diffuse, atmless, direct_axes, diffuse_axes = MakeFluxLookupTable(savename, slope, aspect)

slope    = np.deg2rad(70.0)
aspect   = np.deg2rad(225.0)
for i,savename in enumerate(vnames):        
    direct, diffuse, atmless, direct_axes, diffuse_axes = MakeFluxLookupTable(savename, slope, aspect)

slope    = np.deg2rad(00.0)
aspect   = np.deg2rad(225.0)
for i,savename in enumerate(vnames):        
    direct, diffuse, atmless, direct_axes, diffuse_axes = MakeFluxLookupTable(savename, slope, aspect)
