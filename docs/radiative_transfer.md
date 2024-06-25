# Radiative Transfer

Radiative transfer is key to this project as the solar zenith angles and surface slopes are both high. Incidence angles on the slopes can therefore be close to zero, while the atmospheric path length can be extremely long (>20 times the value from looking vertically for a typical atmospheric scale height on a planet the size of Mars). We use DISORT 4.0.99 to perform radiative transfer calculations accessed through the pyRT.DISORT python library of Kyle Connour[^pyrt]. We operate DISORT in pseudospherical mode due to the prevelance of high solar zenith angles and do the visible and thermal calculations seperately. 

[^pyrt]: The github repository for this library can be found at https://github.com/kconnour/pyRT_DISORT/

The atmophere is described by several parameters that vary with season that can be grouped into:
* Pressure, temperature, dust mixing ratio, and ice mixing ratio as a function of height
* Column-integrated optical depths of dust and ice.

We retrieve these quantities from the Mars Climate Database (MCD)[^mcd] every 15° of Ls and (in the case of the first group above) at 35 vertical levels up to heights of 50km. This retrieval is done for a specific geographic location and martian year.

[^mcd]: Version 6.1 of the Mars Climate database can be found at https://www-mars.lmd.jussieu.fr/ 

Optical properties of dust and water Ice Aerosols come from FITS files provided by Mike Wolff.  These files describe the extinction and scattering cross-sections (single scattering albedo is available from the ratio of these) as well as the phase function moments as a function of particle size and wavelength. Functions provided with pyrt_DISORT combine these optical properties with the MCD outputs of optical depths, vertical distributions of aerosol mixing ratios, and a choice of particle size (we use reff of 1.5 and 2 microns for dust and water ice). Optical depth, single scattering albedo, and particle phase function versus height are calculated and passed to DISORT.

DISORT simulations are done for many (several dozen) wavelengths over the solar or martian blackbody curve (depending on whether we're doing the visible or thermal simulation) and output radiances/fluxes from each separate wavelength co-added. DISORT output includes the direct beam irradiance (in the visible case) as well as diffuse radiance from directions provided as DISORT input. We obtain diffuse DISORT radiances at 360 azimuth angles and 180 zenith angles (from zenith to nadir). Note radiances from below horizonal are from surface scattering/emission, which is invisible to flat surfaces, but will be important for sloping surfaces. These results are generated for many solar zenith angles (or surface temperatures), and several albedo (or emissivity) values. Although all solar azimuths and zenith angles are simulated, the results are only accurate for the Ls of the MCD results used. We loop though these seasons (every 15° of Ls) to update the atmospheric description and adjust the solar flux by the changing Mars-Sun distance. 

Additional python scripts take knowledge of the surface slope and aspect and turn the all-sky radiance results into total diffuse flux onto that surface.  A flat surface will only see radiance scattered downward from the sky, but slopes will also see some reflected/emitted radiation from the surrounding surface. The direct beam flux onto a sloping surface is calculated by renormalizing the DISORT direct bean flux by the ratio of cosine of the incidence angles on flat and sloping terrain i.e. the direct beam passes through the same atmosphere for both flat terrain or slopes, but intersects those two surfaces differently. Direct and diffuse fluxes are combined and reprocessed into two lookup tables:
1. Total visible flux for the surface of interest as a function of season, solar azimuth, zenith angle and surface albedo.
2. Total thermal flux for the surface of interest as a funcion of season, surface temperature, and emissivity.

Note, when using these tables it is the albedo, emissivity, and temperature of the flat terrain that should be used, even when figuring out the flux incident onto the sloping facet.  The surrounding infinite flat plain affects the atmospheric radiation and any sloping facet; however, the sloping facet is small by comparison and does not affect the atmosphere or contribute significant flux to the flat terrain.

## Combining MCD, DISORT, and the IDL thermal model
The IDL thermal model provides the solar azimuth/zeith angle and, for the flat surface, the temperature and albedo/emissivity (which changes seasonally as CO2 frost comes and goes).  Using these values we can interpolate from the two lookup tables described at the end of the Raditive transfer section to get the flux incident on the sloping facet.  

Getting this thermal solution for the flat facet is problematic as we don't know the temperatures and albedo/emissiivty apriori and those quantities affect the radiation incident on the flat terrain. We iterate to converge on the correct answwer.  Fluxes on the flat terrain are first calculated with some simple atmospheric assumuptions, the resulting temperatures and frost behavior are used to update the thermal/visible fluxes and the temperature simulation run again etc... Typically itterating 4-5 times is enough to have temperatures repeat to within a fraction of a degree at all times of year.

Once a self-consistent thermal solution for the flat terrain is known, then that information can be used in the radiative transfer lookup tables to find the fluxes on the sloping surface. No iteration is required on this step as, although the flat surface behavior affects the slope, the slopes temperatures has not effect on the surrounding flat surface.


