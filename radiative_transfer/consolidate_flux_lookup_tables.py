import numpy as np
import pickle


#
# Setup Ls values and all the filenames to consolidate into four tables
#
lsubs              = np.arange(0, 361, 30)
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


# for i,loadname in enumerate(vis_fname):        
#     with open(loadname, 'rb') as f:
#         data = pickle.load(f)
#     if isinstance(data, dict):
#         for key in data:
#             globals()[key] = data[key]
#     del f,key,data

#     if i==0:
#         NAZ, NINC, NALB = diffuse.shape
#         vis = np.zeros((NLS, NAZ, NINC, NALB))
#         vis_axes = (lsubs, solarAZ_i, UMU0_i, ALB_i)

#     vis[i,:,:,:] = diffuse
#     for j in range(NALB):
#         vis[i,:,:,j] = vis[i,:,:,j] + direct


# for i,loadname in enumerate(flat_vis_fname):        
#     with open(loadname, 'rb') as f:
#         data = pickle.load(f)
#     if isinstance(data, dict):
#         for key in data:
#             globals()[key] = data[key]
#     del f,key,data

#     if i==0:
#         NAZ, NINC, NALB = diffuse.shape
#         flat_vis = np.zeros((NLS, NAZ, NINC, NALB))
#         flat_vis_axes = (lsubs, solarAZ_i, UMU0_i, ALB_i)

#     vis[i,:,:,:] = diffuse
#     for j in range(NALB):
#         flat_vis[i,:,:,j] = vis[i,:,:,j] + direct


sname = thermal_fname[0].replace('THERMAL_', '').replace('_LSUBS_' + f"{lsubs[0]:03}",'')
with open(sname, 'wb') as f:
    pickle.dump({
        'thermal': thermal,
        'thermal_axes': thermal_axes,
        'flat_thermal': flat_thermal,
        'flat_thermal_axes': flat_thermal_axes
        # 'vis': vis,
        # 'vis_axes': vis_axes,
        # 'flat_vis': flat_vis,
        # 'flat_vis_axes': flat_vis_axes
    }, f)
