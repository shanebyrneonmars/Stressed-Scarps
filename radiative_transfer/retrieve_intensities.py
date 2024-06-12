import numpy as np
import pyrt
import disort
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle
import time
from scipy.interpolate import interpn
import pandas as pd

start_time = time.time()

rootfolder = '/Users/shane/Desktop/Avalanche/Stressed-Scarps/'


PLANK        = False               # False to do the scattering solution only, no emission. True to include emission
MCD          = False               # True to use MCD output, False to use a simple exponential atmosphere
lsubs        = 30.0
mcd_fname    = rootfolder+'datafiles/MCD_Output/MCD_MY33_83.8N_235E_PTDI.txt'


if PLANK == True:
    NCOL  = 35
    NINC  = 1
    ALB_i = np.array([0.1, 0.05, 0.0])
elif PLANK == False:
    NCOL  = 3 #78
    NINC  = 3 #100
    ALB_i = np.array([0.15,0.25,0.65])


NUMU  = 180
NAZ   = 360
NALB  = ALB_i.shape[0]

# Solar azimuth and cosine of the incidence angle
UMU0_i = np.linspace(1.0,0.003,NINC)    # array of UMU0 to cycle through
UMU0   = 0.707                          # default value of UMU0 to use (only happens when PLANK==True) 
PHI0   = 180.0                          # Solar azimuth angle in degrees

# Return intensities from many user-defined directions in UMU/PHI
ONLYFL = False
USRANG = True
PHI  = np.linspace(0.0+0.5*360.0/NAZ, 360.0-0.5*360.0/NAZ, NAZ)     # Azimuthal angles of results
UMU  = np.linspace(-1.0+0.5*2.0/NUMU, 1.0-0.5*2.0/NUMU, NUMU)       # COS(ELEV) angles of results
UMU  = np.where(np.abs(UMU) < 0.000001, 0.000001, UMU)              # Set minimum abs(UMU) to 0.000001 to avoid numerical issues

# Radiative boundary conditions
FBEAM = 0.0                 # The beam irradiance. This is set in the loop below.
FISOT = 0.0                 # Top Boundary diffuse downward flux - should be zero for an atmosphere
TTEMP = 0.0                 # Top Boundary temperature. Not the temperatue of the uppermost atmospheric layer (best to set to zero for an atmophere)
TEMIS = 1.0                 # Top Boundary emissivity.  Value doesn't matter if TTEMP is zero.
BTEMP = 150.0               # Surface temperature. NOT the temperature of the lowermost atmospheric layer. This will be overwritten later in the PLANK==True case. It's a harmless placeholder in the PLANK==False case
ALBEDO = 0.25               # This is reset in the loop over albedos below. Also 1.0-emissivity of the surface when PLANK==True
LAMBER = True               # Use a Lambert surface - maybe make this fancier someday, but not today

BTEMP_i = np.linspace(140, 300, 161)               # Surface temperatures to cycle through if PLANK is True
NTMP    = BTEMP_i.shape[0]                         # Number of temperatures to loop through

# Setup wavelength arrays
if PLANK == False:
    wvlo = 0.1             # starting wavelength in microns
    wvhi = 4.0             # ending wavelength in microns
elif PLANK == True:
    wvlo = 5.0             # starting wavelength in microns
    wvhi = 40.0            # ending wavelength in microns 

spectral_width = (wvhi - wvlo)/NCOL   # Spectral resolution in microns
wvlen = np.linspace(wvlo + 0.5*spectral_width, wvhi-0.5*spectral_width, NCOL) # wavelengths in microns
WVNMHI = pyrt.wavenumber(wvlen - spectral_width*0.5)         # High wavenumber of each wavelength bin
WVNMLO = pyrt.wavenumber(wvlen + spectral_width*0.5)         # Low wavenumber of each wavelength bin


if MCD == False:
    #Set up atmospheric pressure, temperature, density profiles
    SCALEHEIGHT = 10.8           # Scale height in km
    atm_n = 14                   # Number of atmospheric layers   
    atm_r = 1.1                  # Ratio of geometric progression for the atmospheric layers
    MAXHEIGHT = 7.0*SCALEHEIGHT  # Top of Atmosphere in Scale Heights
    H_LYR     = np.flip(np.insert(np.cumsum(MAXHEIGHT * (1.0 - atm_r) / (1.0 - atm_r**atm_n) * atm_r**(np.arange(0,atm_n,1))),0,0)) # Atmospheric Layer Boundaries (Num of Layers +1)

    TEMPER              = np.linspace(180, 150, num=H_LYR.shape[0])   # create a temperature profile with 14 layers from 150 to 250 K
    pressure_profile    = 610.0 * np.exp(-H_LYR / SCALEHEIGHT)
    column_density      = pyrt.column_density(pressure_profile, TEMPER, H_LYR)

    # Set up Conrath profiles for dust and ice
    altitude_midpoint   = (H_LYR[:-1] + H_LYR[1:]) / 2.0              # mid point of layers
    dust_profile = pyrt.conrath(altitude_midpoint, 1.0, SCALEHEIGHT, 0.1) # Conrath profile normalized to mixing ratio of 1.0 at the surface and nu value of 0.1. Note absolute values shouldn't matter as total optical depth is specified later.
    ice_profile  = pyrt.conrath(altitude_midpoint, 1.0, SCALEHEIGHT, 0.1) # Conrath profile normalized to mixing ratio of 1.0 at the surface and nu value of 0.1. Note absolute values shouldn't matter as total optical depth is specified later.

    NADIRDUSTTAU = 0.1                # Dust Optical Depth - note this should be the EXTINCTION optical depth (not absorption)
    NADIRICETAU  = 0.0                # Ice Optical Depth - note this should be the EXTINCTION optical depth (not absorption)

    if PLANK == True:
        savename = rootfolder+'pickle/DISORT_THERMAL_SimpleAtmosphere_TAUD_'+"{:0.2f}".format(NADIRDUSTTAU)+'_TAUI_'+"{:0.2f}".format(NADIRICETAU)+'.pkl'
    elif PLANK == False:
        savename = rootfolder+'pickle/DISORT_VIS_SimpleAtmosphere_TAUD_'+"{:0.2f}".format(NADIRDUSTTAU)+'_TAUI_'+"{:0.2f}".format(NADIRICETAU)+'.pkl'
elif MCD == True:
    H_LYR, pressure_profile, TEMPER, column_density, dust_profile, ice_profile, NADIRDUSTTAU, NADIRICETAU = process_mcd_output_for_DISORT(lsubs,mcd_fname)

    if PLANK == True:
        savename = rootfolder+'pickle/DISORT_THERMAL_'+mcd_fname.replace('PTDI.txt','LSUBS_').split('/')[-1]+"{:03d}".format(round(lsubs))+'.pkl'
    elif PLANK == False:
        savename = rootfolder+'pickle/DISORT_VIS_'+mcd_fname.replace('PTDI.txt','LSUBS_').split('/')[-1]+"{:03d}".format(round(lsubs))+'.pkl'


NADIRICETAU  = 0.00                # Enforce a low Ice Optical Depth at all times until MCD results are better understood

# Set up particle sizes for dust and ice
dust_particle_size = np.ones(dust_profile.shape[0]) * 1.5
ice_particle_size = np.ones(ice_profile.shape[0]) * 2.0

# Load the scattering data for dust and ice
dwolff = fits.open(rootfolder+'datafiles/Wolff_Aerosols/mars045i_all_area_s0780.fits')
particle_size_grid       = dwolff[5].data     # 25 particle sizes
wavelength_grid          = dwolff[6].data     # 228 wavelength values
pmom_grid                = dwolff[2].data     # 160 scattering moments at each particle size and wavelength
extinction_cross_section = (dwolff[1].data)[:,:,0]
scattering_cross_section = (dwolff[1].data)[:,:,1]
ext_grid                 = pyrt.extinction_ratio(extinction_cross_section, particle_size_grid, wavelength_grid, 9.3) # calculate extinction factor ratio at each wavelength relative to 9.3 microns

dust_ssa           = pyrt.regrid(scattering_cross_section / extinction_cross_section, particle_size_grid, wavelength_grid, dust_particle_size, wvlen)
dust_legendre      = pyrt.regrid(pmom_grid,                                           particle_size_grid, wavelength_grid, dust_particle_size, wvlen)
dust_ext           = pyrt.regrid(ext_grid,                                            particle_size_grid, wavelength_grid, dust_particle_size, wvlen)
dust_optical_depth = pyrt.optical_depth(dust_profile, column_density, dust_ext, NADIRDUSTTAU)

iwolff = fits.open(rootfolder+'datafiles/Wolff_Aerosols/droxtal_050_tmat1_reff_v010_ver121.fits')
particle_size_grid       = iwolff[5].data     # 54 particle sizes
wavelength_grid          = iwolff[6].data     # 445 wavelength values
pmom_grid                = iwolff[2].data     # 192 scattering moments at each particle size and wavelength
extinction_cross_section = (iwolff[1].data)[:,:,0]
scattering_cross_section = (iwolff[1].data)[:,:,1]
ext_grid                 = pyrt.extinction_ratio(extinction_cross_section, particle_size_grid, wavelength_grid, 9.3) # calculate extinction factor ratio at each wavelength relative to 9.3 microns

ice_ssa           = pyrt.regrid(scattering_cross_section / extinction_cross_section, particle_size_grid, wavelength_grid, ice_particle_size, wvlen)
ice_legendre      = pyrt.regrid(pmom_grid,                                           particle_size_grid, wavelength_grid, ice_particle_size, wvlen)
ice_ext           = pyrt.regrid(ext_grid,                                            particle_size_grid, wavelength_grid, ice_particle_size, wvlen)
ice_optical_depth = pyrt.optical_depth(ice_profile, column_density, ice_ext, NADIRICETAU)

rayleigh_co2 = pyrt.rayleigh_co2(column_density, wvlen)
dust_column  = pyrt.Column(dust_optical_depth, dust_ssa, dust_legendre)
ice_column   = pyrt.Column(ice_optical_depth, ice_ssa, ice_legendre)
model        = rayleigh_co2 + dust_column + ice_column

SSALB = model.single_scattering_albedo
PMOM  = model.legendre_coefficients
DTAUC = model.optical_depth
DTAUC = np.where(DTAUC < 0.001, 0.001, DTAUC)  # Set minimum optical depth to 0.00001 to avoid numerical issues

#
# Miscellaneous other variables for DISORT
#

MAXCLY = len(H_LYR)-1
MAXMOM = PMOM.shape[0] - 1
MAXCMU = 32                     # Number of streams. Mike Wolff recommends 16, tests show that it's no different than 64
MAXPHI = len(PHI)
MAXUMU = len(UMU)
MAXULV = 2                      # Two atmospheric levels are set to the top of the atmosphere and the base

ACCUR = 0.0                    # 
DELTAMPLUS = True              # yes to delta-M scaling, although docs say this is on all the time now
DO_PSEUDO_SPHERE = True        # Yes to pseudo-spherical correction
EARTH_RADIUS     = 3376200.0   # Planetarty Radius in the same units as H_LYR - needed for Pseudo-spherical correction
HEADER = ''                    # Info to put in the header - '' should make the header disappear, but it doesn't
PRNT   = [False, False, False, False, False] # Display output of different sorts
IBCND  = False                 # Do not use - this seems to be for overall albedo and transmission calculations as a function of incidence angle

# Two atmospheric levels are set to the top of the atmosphere and the base
USRTAU = True
MAXTAU = np.sum(DTAUC,axis=0)  # Maximum optical depth (at the surface) at each wavelength
UTAU   = [0.0,MAXTAU[0]]       # This will be updated at every wavelength step in the loop below as maximum optical depth changes

# Empty arrays to store the output of DISORT - not actually used in the end
ALBMED = pyrt.empty_albedo_medium(MAXUMU)
FLUP   = pyrt.empty_diffuse_up_flux(MAXULV)
RFLDN  = pyrt.empty_diffuse_down_flux(MAXULV)
RFLDIR = pyrt.empty_direct_beam_flux(MAXULV)
DFDT   = pyrt.empty_flux_divergence(MAXULV)
UU     = pyrt.empty_intensity(MAXUMU, MAXULV, MAXPHI)
UAVG   = pyrt.empty_mean_intensity(MAXULV)
TRNMED = pyrt.empty_transmissivity_medium(MAXUMU)

# More complicated Hapke parameters can be used by setting the following variables
RHOQ  = pyrt.make_empty_rhoq(MAXCMU)
RHOU  = pyrt.make_empty_rhou(MAXCMU)
EMUST = pyrt.make_empty_emust(MAXUMU)
BEMST = pyrt.make_empty_bemst(MAXCMU)
RHO_ACCURATE = pyrt.make_empty_rho_accurate(MAXUMU,MAXPHI)

# Create flux arrays for each wavelength to use as FBEAM later
flux  = np.zeros(NCOL)
if PLANK == False:
    kb    = 1.380649e-23  # Boltzmann's constant
    c     = 2.99792458e8  # Speed of light
    h     = 6.626070e-34  # Planck's constant
    Teff  = 5770.0        # Solar effective temperature
    ww    = wvlen*1e-6    # Wavelength in meters
    rsun  = 6.96e8        # Solar radius in meters
    rau   = 1.496e11      # 1 AU in meters
    aaa   = 1.52366231    # Orbital semi-major axis
    ecc   = 0.09341233    # Orbital eccentricity
    lsp   = 250.870       # Ls of perihelion
 
    sol_dist = rau * aaa * (1.0-ecc**2) / (1.0 + ecc*np.cos(np.deg2rad(lsubs-lsp))) # Solar distance in meters
    bb       = (2.0*h*c**2/(ww**5)) * (1.0/( np.exp(h*c/(ww*kb*Teff))-1.0 ))        # Radiances from the Planck function
    flux     = np.pi * bb * spectral_width*1e-6 * (rsun/sol_dist)**2                # Flux at each wavelength in W/m^2
totalvisflux = np.sum(flux)

# Arrays to store the output as a function of direction, wavelength, and solar incidence angle
uu_surf_inc = np.zeros((NUMU,NAZ,NINC,NALB))
rfldn_inc   = np.zeros((NINC,NALB))
flup_inc    = np.zeros((NINC,NALB))
rfldir_inc  = np.zeros((NINC,NALB))

# Setup two output possibilities for VIS or THERMAL cases. We use the same array structuresm, but VIS uses many INC angles (NINC) and THERMAL uses many BTEMP (NTMP) values 
if PLANK == False:
    uu_surf_inc = np.zeros((NUMU,NAZ,NINC,NALB))
    rfldn_inc   = np.zeros((NINC,NALB))
    flup_inc    = np.zeros((NINC,NALB))
    rfldir_inc  = np.zeros((NINC,NALB))
    j_i         = UMU0_i
elif PLANK == True:
    uu_surf_inc = np.zeros((NUMU,NAZ,NTMP,NALB))
    rfldn_inc   = np.zeros((NTMP,NALB))
    flup_inc    = np.zeros((NTMP,NALB))
    rfldir_inc  = np.zeros((NTMP,NALB))
    j_i         = BTEMP_i

# Call DISORT for every wavelength, solar incidence angle, and albedo
for k,ALBEDO in enumerate(ALB_i):        # Albedo (or emissivity) loop
    for j,jj in enumerate(j_i):          # Loop for over incidence angle (UMU0_i when PLANK==False) or Surface Temperature (BTEMP_i when PLANK==True)
        if PLANK == False:
            UMU0 = jj
        elif PLANK == True:
            BTEMP = jj

        for i,wv in enumerate(wvlen):   # Loop over wavelength
            FBEAM = flux[i]                      # Set the beam irradiance to the flux at the current wavelength (zeros when PLANK==True)
            UTAU  = np.array([0.0,MAXTAU[i]])    # Return radiances at the top and bottom of the atmosphere (bottom TAU is wavelength dependant)

            HEADER = 'Albedo: ' + str(ALBEDO) + ' UMU0: ' + str(UMU0) + ' FBEAM: ' + str(FBEAM)+ ' UTAU: ' + str(UTAU[1])
            rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = \
                disort.disort(USRANG, USRTAU, IBCND, ONLYFL, PRNT, PLANK, LAMBER,
                              DELTAMPLUS, DO_PSEUDO_SPHERE, DTAUC[:, i], SSALB[:, i],
                              PMOM[:, :, i], TEMPER, WVNMLO[i], WVNMHI[i],
                              UTAU, UMU0, PHI0, UMU, PHI, FBEAM, FISOT,
                              ALBEDO, BTEMP, TTEMP, TEMIS, EARTH_RADIUS, H_LYR, RHOQ, RHOU,
                              RHO_ACCURATE, BEMST, EMUST, ACCUR, HEADER, RFLDIR,
                              RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED)
            
            uu_surf_inc[:,:,j,k] = uu_surf_inc[:,:,j,k] + uu[:,1,:]
            rfldir_inc[j,k]      = rfldir_inc[j,k]      + rfldir[1]
            rfldn_inc[j,k]       = rfldn_inc[j,k]       + rfldn[1]
            flup_inc[j,k]        = flup_inc[j,k]        + flup[1]



with open(savename, 'wb') as f:
    pickle.dump({
        'uu_surf_inc': uu_surf_inc,
        'rfldir_inc': rfldir_inc,
        'rfldn_inc': rfldn_inc,
        'flup_inc': flup_inc,
        'PHI': PHI,
        'UMU': UMU,
        'UMU0_i': UMU0_i,
        'PHI0': PHI0,
        'ALB_i': ALB_i,
        'BTEMP_i': BTEMP_i,
        'PLANK': PLANK,
        'totalvisflux': totalvisflux
    }, f)


end_time = time.time()
print("Run time: ", np.round(end_time - start_time,1), "seconds")
