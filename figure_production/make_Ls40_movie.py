
import numpy as np
import pickle
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

savename = 'DISORT_PS_VIS_TAUD_0.05_TAUI_0.00.pkl'
slope    = np.deg2rad(70.0)
aspect   = np.deg2rad(225.0)


#
# Restore data 
#
with open(savename, 'rb') as f:
    data = pickle.load(f)
if isinstance(data, dict):
    for key in data:
        globals()[key] = data[key]
del f,key,data

NUMU = uu_surf_inc.shape[0]
NAZ  = uu_surf_inc.shape[1]
NINC = uu_surf_inc.shape[2]
NALB = uu_surf_inc.shape[3]
NSTP = 120

#
# Get solar coordinates
#
solarAZ, solarINC = MartianSolarDay(NSTP,40.0,84.0)
delta_solarAZ   = np.round((solarAZ*180.0/np.pi - PHI0) / (360.0/NAZ)).astype(int)        # Nearest neighbor azimuthal angle
albedo_arg      = 2


#
# Calculate sky-pixel solid angles and incidence angles
#
UMU_up = 0.5*(np.roll(UMU, 1,axis=0) + UMU)
UMU_dn = 0.5*(np.roll(UMU,-1,axis=0) + UMU)
UMU_up[0]  = -1.0
UMU_dn[-1] = 1.0
inc_flat = np.zeros([NUMU,NAZ])
inc_up   = np.zeros([NUMU,NAZ])
inc_dn   = np.zeros([NUMU,NAZ])
for i in range(NAZ):
    inc_flat[:,i] = np.arccos(-UMU)
    inc_up[:,i]   = np.arccos(-UMU_up)
    inc_dn[:,i]   = np.arccos(-UMU_dn)

phi_flat = np.zeros([NUMU,NAZ])
for i in range(NUMU):
    phi_flat[i,:] = np.deg2rad(PHI)

dphi = 2.0*np.pi/NAZ
dinc = np.abs(inc_up-inc_dn)
sa   = dphi*dinc*0.5*(np.sin(inc_up) + np.sin(inc_dn))

cosi = np.cos(slope)*np.cos(inc_flat) + np.sin(slope)*np.sin(inc_flat)*np.cos(phi_flat-aspect)
cosi = np.where(cosi < 0.0, 0.0, cosi)

cosisa = cosi*sa

#
# Loop through 
#

k=0
uu_AZ_corrected = np.roll(uu_surf_inc, delta_solarAZ[k], axis=1)                                        # Roll the sky radiance array in azimuth to match the solar azimuth
solar_INC_arg   = np.abs(np.cos(solarINC[k])-UMU0_i).argmin()
sky_radiance    = uu_AZ_corrected[:,:,solar_INC_arg,albedo_arg]
sky_flux        = np.ma.masked_where(cosi == 0, sky_radiance*cosisa)



#
# Create an animation of the solar power incident on a surface with horizon masking
#

# Create a custom colormap
cmap = plt.cm.turbo
points = []
cmap.set_bad('gray')  # set color for zero values

# Mask zero values
k=0
uu_AZ_corrected = np.roll(uu_surf_inc, delta_solarAZ[k], axis=1)                                        # Roll the sky radiance array in azimuth to match the solar azimuth
solar_INC_arg   = np.abs(np.cos(solarINC[k])-UMU0_i).argmin()
sky_radiance    = uu_AZ_corrected[:,:,solar_INC_arg,albedo_arg]
sky_flux        = np.ma.masked_where(cosi == 0, sky_radiance*cosisa)

data = sky_flux
mname = 'Flux from each pixel (1x1 degrees) [W m$^{-2}$]'
lims = [0.001,1.0]

data  = sky_radiance
mname = 'Radiance [W m$^{-2}$ sr$^{-1}$]'
lims  = [1,10000]


# Make plot with first frame
fig, ax = plt.subplots(figsize=(5, 4))
norm = colors.LogNorm(vmin=lims[0], vmax=lims[1])
im = ax.imshow(data, extent=[0, 360, -1, 1], aspect=90, norm=norm, cmap=cmap)
ax.set_xlabel('Azimuth')
ax.set_ylabel('Sine of Elevation')
ax.xaxis.set_ticks(np.linspace(0,360,9))
plt.tight_layout()

# Add solar spot, keep track of the points for later removal
point = ax.plot(solarAZ[k]*180.0/np.pi,np.cos(solarINC[k]), 'ro',markersize=2)
points.append(point[0])

# Add colorbar
colorbar = fig.colorbar(im, ax=ax, location='top', orientation='horizontal')
colorbar.set_label(mname)

# define function that updates the plot
def update(k):
    if points:
        points[0].remove()
        points.pop(0)
    uu_AZ_corrected = np.roll(uu_surf_inc, delta_solarAZ[k], axis=1)                                        # Roll the sky radiance array in azimuth to match the solar azimuth
    solar_INC_arg   = np.abs(np.cos(solarINC[k])-UMU0_i).argmin()
    sky_radiance    = uu_AZ_corrected[:,:,solar_INC_arg,albedo_arg]
    sky_flux        = np.ma.masked_where(cosi == 0, sky_radiance*cosisa)
    
    data = sky_flux
    data = sky_radiance
    im.set_array(data)

    point = ax.plot(solarAZ[k]*180.0/np.pi,np.cos(solarINC[k]), 'ro',markersize=2)
    points.append(point[0])

    return im,

# Define and save the animation
ani = animation.FuncAnimation(fig, update, frames=range(NSTP), blit=True)
ani.save('/Users/shane/Desktop/movie_rad.mp4', writer='ffmpeg')
plt.show()



