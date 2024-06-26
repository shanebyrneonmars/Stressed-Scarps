# Avalanche Scarp Stress Calculation Outline
Thermoelastic stress in these scarps can cause cracks and sheeting joints that result in exfoliation of slab-like ice blocks. Simulations of this stress can be broken into three main sections:
1. Radiative transfer of the martian atmosphere to determine the energy delivered to the surface.
2. Diffusion of thermal energy in the subsurface to determine how the surface temperature evolves.
3. Thermal expansion and contraction of the ice combined with viscous processes that leads to surface-parallel stress.

The general outline of what happens at each step is described here. The [workflow](Workflow.md "Workflow") document contains details of how to implement it.

## Radiative Transfer
Radiative transfer is key to this project as the solar zenith angles and surface slopes are both high. Incidence angles on the slopes can therefore be close to zero, while the atmospheric path length can be extremely long (>20 times the value from looking vertically for a typical atmospheric scale height on a planet the size of Mars). We use DISORT 4.0.99 to perform radiative transfer calculations accessed through the pyRT.DISORT python library of Kyle Connour[^pyrt]. We operate DISORT in pseudospherical mode due to the prevalence of high solar zenith angles and do the visible and thermal calculations separately. 

[^pyrt]: The github repository for this library can be found at https://github.com/kconnour/pyRT_DISORT/

The atmophere is described by several parameters that vary with season that can be grouped into:
* Pressure, temperature, dust mixing ratio, and ice mixing ratio as a function of height
* Column-integrated optical depths of dust and ice.

We retrieve these quantities from the Mars Climate Database (MCD)[^mcd] every 15° of Ls and (in the case of the first group above) at 35 vertical levels up to heights of 50km. This retrieval is done for a specific geographic location and martian year. Dust optical depths are returned as extinction (not absorption) optical depths specific to the surface pressure at that point in the visible (we divide by 2 to convert to the value at 9.3 &mu;m). Water ice is returned as kg m<sup>-2</sup>, which we convert to optical depth at 9.3 &mu;m by assuming non-overlapping particles and a Q<sub>ext</sub> of 2.62.

[^mcd]: Version 6.1 of the Mars Climate database can be found at https://www-mars.lmd.jussieu.fr/ 

Optical properties of dust and water Ice Aerosols come from FITS files provided by Mike Wolff[^wolff_aerosols].  These files describe the extinction and scattering cross-sections (single scattering albedo is available from the ratio of these) as well as the phase function moments as a function of particle size and wavelength. Functions provided with pyrt_DISORT combine these optical properties with the MCD outputs of optical depths, vertical distributions of aerosol mixing ratios, and a choice of particle size (we use r<sub>eff</sub> of 1.5 and 2 microns for dust and water ice). Optical depth, single scattering albedo, and particle phase function versus height are calculated and passed to DISORT.

[^wolff_aerosols]: Available at https://gemelli.spacescience.org/~wolff/files/aerosols/

DISORT simulations are done for many (several dozen) wavelengths over the solar or martian blackbody curve (depending on whether we're doing the visible or thermal simulation) and output radiances/fluxes from each separate wavelength co-added. DISORT output includes the direct beam irradiance (in the visible case) as well as diffuse radiance from directions provided as DISORT input. We obtain diffuse DISORT radiances at 360 azimuth angles and 180 zenith angles (from zenith to nadir). Note radiances from below horizontal are from surface scattering/emission, which is invisible to flat surfaces, but will be important for sloping surfaces. These results are generated for many solar zenith angles (or surface temperatures), and several albedo (or emissivity) values. Although all solar azimuths and zenith angles are simulated, the results are only accurate for the Ls of the MCD results used. We loop though these seasons (every 15° of Ls) to update the atmospheric description and adjust the solar flux by the changing Mars-Sun distance. 

Additional python scripts take knowledge of the surface slope and aspect and turn the all-sky radiance results into total diffuse flux onto that surface.  A flat surface will only see radiance scattered downward from the sky, but slopes will also see some reflected/emitted radiation from the surrounding surface. The direct beam flux onto a sloping surface is calculated by renormalizing the DISORT direct bean flux by the ratio of cosine of the incidence angles on flat and sloping terrain i.e. the direct beam passes through the same atmosphere for both flat terrain or slopes, but intersects those two surfaces differently. Direct and diffuse fluxes are combined and reprocessed into two lookup tables:
1. Total visible flux for the surface of interest as a function of season, solar azimuth, zenith angle and surface albedo.
2. Total thermal flux for the surface of interest as a function of season, surface temperature, and emissivity.

Note, when using these tables it is the albedo, emissivity, and temperature of the flat terrain that should be used, even when figuring out the flux incident onto the sloping facet.  The surrounding infinite flat plain affects the atmospheric radiation and any sloping facet; however, the sloping facet is small by comparison and does not affect the atmosphere or contribute significant flux to the flat terrain.

## Thermal model description

This project utilizes a 1D thermal diffusion model to simulate temperature at the surface and at depth. The model supports multiple solar system bodies (incl. user defined orbital elements), changes in thermophysical properties with depth, and arbitrary terrain slope and aspect.  In the case of Mars, atmospheric pressure varies seasonally and user-defined elevation is used to calculate the seasonally varying CO<sub>2</sub> frost point.

Slopes receive energy from the sky (both direct and diffusely scattered) and the surrounding terrain. In this model we assume the surrounding terrain is an infinite flat plain that has Lambert scattering properties. As this is a 1D model, sloping surfaces need two model runs. First a model for flat terrain is run, which creates two results files:
1. The results for the flat surface: temperature, CO<sub>2</sub> frost mass etc...
2. The upwelling flux from the flat surface (if requested).
This second output file can then be utilized in the next thermal model run for a sloping surface.  Upwelling reflected visible flux and emitted thermal flux is added to the direct and atmospheric flux 

There are a few steps to running the model for either sloping or flat terrain.  Firstly, the model parameters are chosen and stored in an IDL structure. This is easiest to start by returning a structure with all the necessary fields defined that can be updated using the `return_blank_thermal_model` function with no arguments. Secondly, this structure is passed to a function called `preprocess_thermal_model`, which calculates all the incident radiation and the vertical position of all subsurface layers and their properties. This function saves a file with these quantities (using the IDL save format), which is ready to use as input to the thermal model itself and returns the filename.  Finally, the name of this file is passed to the thermal model function `thermal_model_v07`.  The thermal model is agnostic to the body being simulated and just implements the output of the preprocess_thermal_model function. The thermal model function saves a results file (again with IDL save file format) and returns the filename of the results file. 

If running the thermal model on its own, the assumptions about atmospheric radiation are simplistic and conventional (e.g. scattered radiation is a fixed proportion of the direct beam and atmospheric thermal emission is a fixed proportion of the noon-time direct beam - these proportions are configurable in the parameters structure). The workflow can be summarized in the figure below. To integrate the radiative transfer solutions of DISORT the user intervenes after the preprocessing step and before the thermal model run in a procedure described below.

 ![Thermal model flow chart](thermal_processing_logic.jpg)

## Combining MCD, DISORT, and the IDL thermal model
The IDL thermal model provides the solar azimuth/zenith angle and, for the flat surface, the temperature and albedo/emissivity (which changes seasonally as CO<sub>2</sub> frost comes and goes).  Using these values we can interpolate from the two lookup tables described at the end of the Radiative transfer section to get the flux incident on the sloping facet.  

Getting this thermal solution for the flat facet is problematic as we don't know the temperatures and albedo/emissivity apriori and those quantities affect the radiation incident on the flat terrain. We iterate to converge on the correct answer.  Fluxes on the flat terrain are first calculated with some simple atmospheric assumptions, the resulting temperatures and frost behavior are used to update the thermal/visible fluxes and the temperature simulation run again etc... Typically iterating 4-5 times is enough to have temperatures repeat to within a fraction of a degree at all times of year.

Once a self-consistent thermal solution for the flat terrain is known, then that information can be used in the radiative transfer lookup tables to find the fluxes on the sloping surface. No iteration is required on this step as, although the flat surface behavior affects the slope, the slopes temperatures has no effect on the surrounding flat surface.

## Turning Thermal results into Stresses
We follow the approach of Mellon[^Mellon] to solve for the time varying stress in an initially unfractured viscoelastic solid. No lateral strain can occur, so surface-parallel thermal expansion and contraction is opposed by elastic stresses over short timescales that decay over longer timescales due to grain-size-dependent viscous effects[^Goldsby]. The thermal history at each depth can be used to calculate these stresses, as well as surface-normal displacement. We use a Zenner pinning approach[^Durand] with NPLD dust abundances derived from RADAR sounding data to constrain ice grain sizes to be 10–1000 um.

[^Mellon]: Mellon, M.T., 1997. Small-scale polygonal features on Mars: Seasonal thermal contraction cracks in permafrost. J. Geophys. Res. 102, 25617.
[^Durand]: Durand, G., Weiss, J., Lipenkov, V., Barnola, J.M., Krinner, G., Parrenin, F., Delmonte, B., Ritz, C., Duval, P., Röthlisberger, R., Bigler, M., 2006. Effect of impurities on grain growth in cold ice sheets. J. Geophys. Res. 111, F01015.
[^Goldsby]: Several stress-dependent deformation mechanisms are combined. See: Goldsby, D.L., Kohlstedt, D.L., 2001. Superplastic deformation of ice: Experimental observations. J. Geophys. Res. 106, 11017.

The thermal model above output a variable called `sav_allt` (temperatures at all times and depths) that we will use to calculate expansion and contraction of the ice and the stress this causes. Each depth is independently calculated since the surface-parallel strain is zero. The zero-stress baseline temperature is calculated to be the  annual-average temperature at each depth.
