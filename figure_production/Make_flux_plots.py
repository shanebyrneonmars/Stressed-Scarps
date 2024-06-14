import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interpn

savename = '../pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_030.pkl'
slope    = np.deg2rad(70.0)
aspect   = np.deg2rad(225.0)

direct, diffuse, atmless = MakeFluxLookupTable(savename, slope, aspect)

solarAZ, solarINC = MartianSolarDay(120,40.0,84.0)
cc = np.stack((solarAZ, solarINC))

solarAZ_i       = np.linspace(0,359,360) * np.pi/180.0
r1              = interpn((solarAZ_i,UMU0_i),direct,np.transpose(np.stack((solarAZ, np.cos(solarINC)))))
r2              = interpn((solarAZ_i,UMU0_i),diffuse[:,:,2],np.transpose(np.stack((solarAZ, np.cos(solarINC)))))
r3              = interpn((solarAZ_i,UMU0_i),atmless,np.transpose(np.stack((solarAZ, np.cos(solarINC)))))

hrs = np.linspace(0,24,120) 

fig, ax = plt.subplots()
ax.set_title('Slope of '+str(np.rad2deg(slope))+' Aspect of '+str(np.rad2deg(aspect)))

ax.plot(hrs,r3/totalvisflux, label='No Atmosphere or Surface Scattering',color='black',linestyle='dotted')
ax.plot(hrs,r1/totalvisflux, label='Direct Beam with Adsorption',color='red')
ax.plot(hrs,r2/totalvisflux, label='Atmospheric+Surface Scattering',color='green')
ax.plot(hrs,(r1+r2)/totalvisflux, label='Total',color='blue')

ax.legend()
ax.set_xlim([0,24])
ax.set_ylim([0,1])
ax.set_xlabel('Hour Angle')
ax.set_ylabel('Flux')
ax.xaxis.set_ticks([0,3,6,9,12,15,18,21,24])
plt.tight_layout()
plt.savefig('/Users/shane/Desktop/summary_tau0.jpg', format='jpg')
plt.show()



