
# Avalanche Scarp Stress Calculation Outline

Here we step through the different programs that need to be run to produce some results. There's not much detail of what's happening at each step so readers are encouraged to read the [outline](Outline.md "Outline").

## Radiative Transfer

Several python programs need to be called to generate the flux lookup tables for flat and sloping terrain at each solar position, season, and model atmosphere.

### Retrieve atmospheric parameters from the Mars Climate Database (MCD)

The IDL script `radiative_transfer/MCD_results_retrieval.idl` needs to be run to collect atmospheric information. This is most easily done on the IDL command line as:
```idl
.run radiative_transfer/MCD_results_retrieval.idl
```
The script will use the MCD web interface and the command line program `wget` to retrieve output in a two-step process. The first step retrieves the MCD results page, which contains a link to a text table containing the results. A second `wget` command retrieves this txt table (the original downloaded webpage is removed after striping out the text file link).
We loop over Mars Years 24 to 35 and retrieve two results files from each year. The first provides the pressure, temperature, and dust/ice mixing ratios as a function of season and height, while the second provides dust/ice optical depth/column abundances as a function of season. Filenames generated for MY24 are:
```
MCD_MY24_83.8N_235E_PTDI.txt
MCD_MY24_83.8N_235E_TAU.txt
```
Twelve Mars years with two files each mean 24 text files are expected. These should be stored in `datafiles/MCD_Output/`

The python function `read_mcd_output` can read these text files assuming the tables in each one all have the same format.  This function is contained in the `radiative_transfer/pyRT_DISORT_MCD_utils.py` file and can be imported and run with:
```python
from pyRT_DISORT_MCD_utils import read_mcd_output
tname, cname, rname, cols, rows, vals = read_mcd_output('MCD_MY24_83.8N_235E_PTDI.txt')
```
Table names will be in the array `tname`, column names in `cname`, row names in `rname`, column values (elevation in this case) will be in `cols`, row values (L<sub>s</sub> in this case) will be in `rows`, and the values for all tables in the 3D array `vals`.

We usually need to extract the state of the atmosphere at a specific L<sub>s</sub> and return them in arrays that DISORT can understand. For this we use the python function `process_mcd_output_for_DISORT` in the `radiative_transfer/pyRT_DISORT_MCD_utils.py` file.  It can be run with via:
```python
from pyRT_DISORT_MCD_utils import process_mcd_output_for_DISORT, read_mcd_output
H_LYR, pressure_profile, TEMPER, column_density, dust_profile, ice_profile, NADIRDUSTTAU, NADIRICETAU = process_mcd_output_for_DISORT(lsubs,'MCD_MY24_83.8N_235E_PTDI.txt')
```
Note the filename supplied contains the 'PTDI' string, the function will read this and look for the second file with the 'TAU' string in the same directory. The outputs `H_LYR, pressure_profile, TEMPER, column_density, dust_profile, ice_profile, NADIRDUSTTAU, NADIRICETAU` are respectively the height, pressure, temperature, column density, dust mixing ratio (normalized to the surface value), ice mixing ratio (normalized to the surface value), column dust optical depth, and column ice optical depth (optical depths are extinction, pressure corrected, and at 9.3 &mu;m). These arrays are ordered from the highest elevation to the lowest as DISORT expects. The L<sub>s</sub> in the function call can be any value and interpolation between MCD outputs will be used.

### Running DISORT and saving raw output files
The python program `radiative_transfer/retrieve_intensities.py` cycles through all seasons and calls DISORT to calculate the radiance incident from every direction.  It loops through one atmospheric state file (every 15° of L<sub>s</sub>) that was extracted from the MCD in the previous step (assuming that `MCD == True`; if `MCD == False` then it uses a simple atmosphere - see code for details). It calls `process_mcd_output_for_DISORT` (described above) to get the atmospheric state at each season. The Mars-Sun distance is also updated each season when setting the top of atmosphere flux (`FBEAM`).

DISORT either does the thermal or visible calculation depending on whether `THERMAL == True` or `THERMAL == False`. Both are required so `retrieve_intensities.py` needs to be run twice with different values of this flag. DISORT runs each wavelength separately and the scattering properties of the aerosols vary strongly with wavelength (as does the radiance from the solar and martian Planck curves). `retrieve_intensities.py` runs several dozen wavelengths in both cases and adds the resulting diffuse radiances (and direct beam fluxes in the case of the visible runs). Results are stored in python pickle files for convenience. Each Ls (24 values) produces a results file for both VIS and THERMAL runs (so 48 files produced in total). Filenames have the format:
```
pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_000.pkl
pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_000.pkl
etc...
```
**WARNING:** Running `retrieve_intensities.py` with `THERMAL == False` is very time consuming. On my laptop it takes ~100 minutes for each L<sub>s</sub> value. There are many inefficiencies here that I'll correct in the future - chief among them is that the code simulates all solar incidence angles at each L<sub>s</sub>; however only certain incidence angles are possible for a given L<sub>s</sub> and latitude.

### Preparing DISORT fluxes
We calculated radiances from all directions above, but any given surface sees only half of these directions and some portion of these directions are seen only at grazing angles.  We use the function `MakeFluxLookupTable` contained in the `radiative_transfer/pyRT_DISORT_MCD_utils.py` file to calculate total direct and diffuse flux upon the surface facet of interest and so reduce the dimensionality of the stored results. It can be called with:
```python
direct, diffuse, atmless, direct_axes, diffuse_axes = MakeFluxLookupTable(filename, slope, aspect)
```
The Filename argument is one of the pickle files described above that was output by `retrieve_intensities.py`. The returned arrays, `direct, diffuse, atmless`, contain the direct beam (visible only), diffuse flux (both atmospheric and scattered from the surface), and atmosphere-less direct beam (for comparisons sake).  The `direct_axes` collection contains the solar azimuth and incidence angles i.e. the axes of the 2D arrays contained in `direct` and `atmless`. The `diffuse_axes` collection contains the solar azimuth, incidence angles, and albedo (if working in the visible) or the solar azimuth, surface temperature, and emissivity (if working in the thermal - solar azimuth is only included here to keep the array three-dimensional, it doesn't actually affect the thermal flux).

Each of the 48 output files from `retrieve_intensities.py` needs to be processed in this way and we need results for both flat terrain as well as at least one slope. As well as directly returning into python variables, the `MakeFluxLookupTable` also saves another pickle file with the results (thus generating 96 new pickle files). For one value of L<sub>s</sub>, filenames contain the slope and aspect values and have the format:
```
pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_000_S00_A225.pkl
pickle/DISORT_THERMAL_MCD_MY33_83.8N_235E_LSUBS_000_S70_A225.pkl
pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_000_S00_A225.pkl
pickle/DISORT_VIS_MCD_MY33_83.8N_235E_LSUBS_000_S70_A225.pkl
etc...
```
We loop through all 96 calls with some Python code in the file `radiative_transfer/apply_slopes_to_DISORT_results.py`.

Finally, having results spread over many pickle files doesn't make life easy when using with a thermal model written in IDL. So we run the Python code in `radiative_transfer/consolidate_flux_lookup_tables.py` to collapse all these files for different seasons into a single file and save it as a FITS format file. FITS is a common astronomical format that can store different data structures and has the wonderful property of being readable/writeable by both Python and IDL. VIS and Thermal results for sloping and flat terrain are stored for all seasons in the same file. Thus all our radiative transfer work boils down to arrays stored in a single file e.g.
```
pickle/DISORT_MCD_MY33_83.8N_235E_S70A225.fits
```

## Thermal model 


### Filename conventions
The model parameter structure (`sp`) has a field called (`sp.run_number`) that will be applied to exported files at various steps.  It's useful to give this string an informative name 

## Stress model


## Plots of outputs
