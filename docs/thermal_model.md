# Thermal model description

This project utilizes a 1D thermal diffusion model to simulate temperature at the surface and at depth. The model supports multiple solar system bodies (incl. user definged orbital elements), changes in thermophysical properties with depth, and arbitrary terrain slope and aspect.  In the case of Mars, atmospheric pressure varies seasonally and user-defined elevation is used to calculate the seasonally varying CO<sub>2</sub> frost point.

Slopes receive energy from the sky (both direct and diffusely scattered) and the surrounding terrain. In this model we assume the surrounding terrain is an infinite flat plain that has Lambert scattering properties. As this is a 1D model, sloping surfaces need two model runs. First a model for flat terrain is run, which creates two results files:
1. The results for the flat surface: temperature, CO<sub>2</sub> frost mass etc...
2. The upwelling flux from the flat surface (if requested).
This second output file can then be utilized in the next thermal model run for a sloping surface.  Upwelling reflected visible flux and emitted thermal flux is added to the direct and atmospheric flux 

There are a few steps to running the model for either sloping or flat terrain.  Firstly, the model parameters are chosen and stored in an IDL structure. This is easiest to start by returning a structure with all the necessary fields defined that can be updated using the `return_blank_thermal_model` function with no arguements. Secondly, this structure is passed to a function called `preprocess_thermal_model`, which calculates all the incident radiation and the vertical position of all subsurface layers and their properties. This function saves a file with these quantities (using the IDL save format), which is ready to use as input to the termal model itself and returns the filename.  Finally, the name of this file is passed to the thermal model function `thermal_model_v07`.  The thermal model is agnostic to the body being simulated and just implements the output of the preprocess_thermal_model function. The thermal model function saves a results file (again with IDL save file format) and returns the filename of the results file. 

## Filename conventions
The model parameter structure (`sp`) has a field called (`sp.run_number`) that will be applied to exported files at various steps.  It's useful to give this string an informative name 



 ![Thermal model flow chart](thermal_processing_logic.jpg)

