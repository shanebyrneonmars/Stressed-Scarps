# Thermal model descrition

This project utilizes a 1D thermal diffusion model to simulate temperature at the surface and at depth. The model supports multiple solar system bodies (incl. user definged orbital elements), changes in thermophysical properties with depth, and arbitrary terrain slope and aspect.  In the case of Mars, atmospheric pressure varies seasonally and user-defined elevation is used to calculate the seasonally varying CO ~2~ frost point.

As this is a 1D model, sloping surfaces need two model runs. First a model for flat terrain is run, which creates two results files:
1. The results for the flat surface: temperature, CO2 frost mass etc...
2. The upwelling flux from the flat surface (if requested).
This second output file can then be utilized in the next thermal model run for a sloping surface. Slopes are assumed to be surrounded by infinite flat planes that behave as Lambert surfaces.  Upwelling reflected visible flux and emitted thermal flux is added to the direct and atmospheric flux 

 ![Thermal model flow chart](thermal_processing_logic.jpg)

