# emit-sds-l2b

Welcome to the EMIT Level 2b science data system repository.  To understand how this repository is linked to the rest of the emit-sds repositories, please see [the repository guide](https://github.jpl.nasa.gov/emit-sds/emit-main/wiki/Repository-Guide).

The repository supports two methods for creating mineral maps.  The first method runs and post-processes the output of tetracorder, aggregating the band depth maps from individual library spectra into a combined 10-channel product that shows the band depths of the 10 EMIT mineral classes. This occurs in two main steps: 

1) Generating mineral maps via tetracorder
2) Post-processing these maps to the EMIT-10 estimated spectral abundance

These steps can be run in a coordinated fashion using:

> sh run_tetracorder.sh [tetracorder_output_directory] [aggregated_output_base] [l2a_reflectance_file] [l2a_reflectance_uncertainty_file] 

Before running, compile the program "specpr2envi" and use it to convert the tetracorder "specpr" format library into an equivalent ENVI spectral library.  Extra header fields contain critical ancillary information like the record numbers and names.  This step must be repeated for each reference library used (multiple are supported). Syntax:
> specpr2envi [specpr_library_path] [output_envi_library_path] 


