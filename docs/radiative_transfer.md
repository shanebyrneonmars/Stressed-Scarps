# Radiative Transfer

We use DISORT 4.0.99 to perform radiative transfer calculations accessed through the pyRT.DISORT python library of Kyle Connour. We operate DISORT in pseudospherical mode and do the visible and thermal calculations seperately.

The atmophere is described by several parameters that vary with season:
* Pressure, temperature, dust mixing ratio, and ice mixing ratio as a function of height
* Column-integrated optical depths of dust and ice.

We retrieve these quantities from the Mars Climate Database (MCD) every 30° of Ls and (in the case f the first bullet above) at 35 vertical levels up to heights of 50km. This retrieval is done for a specific geographic location and martian year.

Optical properties of dust and water Ice Aerosols come from FITS files provided by Mike Wolff.  These files describe the extinction and scattering cross-sections (single scattering albedo is available from the ratio of these) as well as the phase function moments as a function of particle size and wavelength. Functions provided with pyrt_DISORT combine these optical properties with  MCD outputed optical depths, vertical distribution of mixing ratio, and a choice of particle size (we use reff of 1.5 and 2 microns for dust and water ice). Optical depth, single scattering albedo, and particle phase function versus height is calculated and passed to DISORT.

DISORT simulations are done for many (several dozen) wavelengths over the solar or martian blackbody curve (depending on whether we're doing the visible or thermal simulation) and output radiances/fluxes from each separate wavelength added. Output includes the direct beam irradiance (in the visible case) as well as radiance from directions provided as DISORT input. We obtain DISORT radiances at 360 azimuth angles and 180 zenith angles (from zenith to nadir). Note radiances from below horizonal are from surface scattering/emission, which is invisible to flat surfaces, but will be important for sloping surfaces. These results are generated for many solar zenith angles (or surface temperatures), and several albedo (or emissivity) values. Although all solar azimuths and zenith angles are simulated, the results are only accurate for the Ls of the MCD results used. We loop though these seasons (every 15° of Ls) to update the atmospheric description and adjust the solar flux by the changing Mars-Sun distance. 

Additional python scripts take knowledge of the surface slope and aspect and turn the all-sky radiance results into total diffuse flux onto the surface.  A flat surface will only see radiance scattered downward from the sky, but slopes will see reflected/emitted radiation from the surrounding surface. The direct beam flux onto a sloping surface is calculated by renormalizing the DISORT results by the ratio of cosine of the incidence angles on flat and sloping terrain i.e. the direct beam passes through the same atmosphere for both flat terrain or slopes, but intersects those two surfaces differently. Output from this stage is reprocessed into two lookup tables:
1. Total visible flux for the surface of interest as a function of season, solar azimuth, zenith angle and surface albedo.
2. Total thermal flux for the surface of interest as a funcion of season, surface temperature, and emissivity.

The IDL thermal model provides the solar azimuth/zeith-angle and, for the flat surface, the temperature and albedo/emissivity (which changes seasonally as CO2 frost comes and goes).  Using these values we can interpolate from the two lookup tables described above the to get the flux incident on the sloping facet.  Flux for the flat facet is more problematic as we don't know the temperatures and albedo/emissiivty apriori. In this case we iterate to converge on the correct answwer.  Fluxes are first calculated with some simple atmospheric assumuptions, the resulting temperatures and frost behavior is used to update the fluxes and the temperature recalculated etc... Typically itterating 4-5 times is enough to have temperatures repeat to within a fraction of a degree at all times of year.


