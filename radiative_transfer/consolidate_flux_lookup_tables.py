import numpy as np
import pickle
from astropy.io import fits

#
# Setup Ls values and all the filenames to consolidate into four tables
#
lsubs              = np.arange(0, 361, 15)
NLS                = lsubs.shape[0]
thermal_fname      = ['../pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_' + f"{item:03}" + '_S70A225.pkl' for item in lsubs]
flat_thermal_fname = ['../pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_' + f"{item:03}" + '_S00A225.pkl' for item in lsubs]
vis_fname          = ['../pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_'     + f"{item:03}" + '_S70A225.pkl' for item in lsubs]
flat_vis_fname     = ['../pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_'     + f"{item:03}" + '_S00A225.pkl' for item in lsubs]

for i,loadname in enumerate(thermal_fname):        
    with open(loadname, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    if i==0:
        NAZ, NTMP, NEMM = diffuse.shape
        thermal = np.zeros((NLS, NTMP, NEMM))
        thermal_axes = (lsubs, BTEMP_i, 1.0-ALB_i)


    thermal[i,:,:] = diffuse[0,:,:]


for i,loadname in enumerate(flat_thermal_fname):        
    with open(loadname, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    if i==0:
        NAZ, NTMP, NEMM = diffuse.shape
        flat_thermal = np.zeros((NLS, NTMP, NEMM))
        flat_thermal_axes = (lsubs, BTEMP_i, 1.0-ALB_i)

    flat_thermal[i,:,:] = diffuse[0,:,:]


for i,loadname in enumerate(vis_fname):        
    with open(loadname, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    if i==0:
        NAZ, NINC, NALB = diffuse.shape
        vis = np.zeros((NLS, NAZ, NINC, NALB))
        vis_axes = (lsubs, solarAZ_i, UMU0_i, ALB_i)

    vis[i,:,:,:] = diffuse
    for j in range(NALB):
        vis[i,:,:,j] = vis[i,:,:,j] + direct


for i,loadname in enumerate(flat_vis_fname):        
    with open(loadname, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    if i==0:
        NAZ, NINC, NALB = diffuse.shape
        flat_vis = np.zeros((NLS, NAZ, NINC, NALB))
        flat_vis_axes = (lsubs, solarAZ_i, UMU0_i, ALB_i)

    flat_vis[i,:,:,:] = diffuse
    for j in range(NALB):
        flat_vis[i,:,:,j] = flat_vis[i,:,:,j] + direct


sname = thermal_fname[0].replace('THERMAL_', '').replace('_LSUBS_' + f"{lsubs[0]:03}",'')
with open(sname, 'wb') as f:
    pickle.dump({
        'thermal': thermal,
        'thermal_axes': thermal_axes,
        'flat_thermal': flat_thermal,
        'flat_thermal_axes': flat_thermal_axes,
         'vis': vis,
         'vis_axes': vis_axes,
         'flat_vis': flat_vis,
         'flat_vis_axes': flat_vis_axes
    }, f)



hdr = fits.Header()
hdr['COMMENT1'] = "Thermal and flat_thermal are 3D arrays with dimensions (Ls, Temperature, Emissivity)"
hdr['COMMENT2'] = "Vis and flat_vis are 4D arrays with dimensions (Ls, solar Azimuth, Solar Incidence, Albedo)"
empty_primary = fits.PrimaryHDU(header=hdr)
image_hdu1 = fits.ImageHDU(data=thermal,      name="Thermal flux on slope")
image_hdu2 = fits.ImageHDU(data=flat_thermal, name="Thermal flux on flat surface")
image_hdu3 = fits.ImageHDU(data=thermal_axes[0].astype(float), name="Ls")
image_hdu4 = fits.ImageHDU(data=thermal_axes[1].astype(float), name="Temperature")
image_hdu5 = fits.ImageHDU(data=thermal_axes[2].astype(float), name="Emmissivity")
image_hdu6 = fits.ImageHDU(data=vis_axes[1].astype(float), name="Solar Azimuth")
image_hdu7 = fits.ImageHDU(data=vis_axes[2].astype(float), name="Cos of Solar Incidence")
image_hdu8 = fits.ImageHDU(data=vis_axes[3].astype(float), name="Albedo")
image_hdu9 = fits.ImageHDU(data=vis,      name="Visible flux on slope")
image_hdu10= fits.ImageHDU(data=flat_vis, name="Visible flux on flat surface")
hdul = fits.HDUList([empty_primary, image_hdu1, image_hdu2, image_hdu3, image_hdu4, image_hdu5, image_hdu6, image_hdu7, image_hdu8, image_hdu9, image_hdu10])
hdul.writeto(sname.replace('.pkl', '.fits'),overwrite=True)
