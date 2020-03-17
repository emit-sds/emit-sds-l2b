# emit-sds-l2b

Welcome to the EMIT Level 2b science data system repository.  To understand how this repository is linked to the rest of the emit-sds repositories, please see [the repository guide](https://github.jpl.nasa.gov/emit-sds/emit-main/wiki/Repository-Guide).

The repository supports two methods for creating mineral maps.  The first method runs and post-processes the output of tetracorder, aggregating the band depth maps from individual library spectra into a combined 10-channel product that shows the band depths of the 10 EMIT mineral classes. This occurs in two main steps: 

1) Generating mineral maps via tetracorder
2) Post-processing these maps to the EMIT-10 estimated spectral abundance

Step 1) can be executed with:

> sh run_tetracorder.sh [l2a_reflectance_file] [scaled_reflectance_destination] [tetracorder_output_directory]

where the latter two arguments are output locations.

Before postprocessing tetracorder output, compile the program "specpr2envi" and use it to convert the tetracorder "specpr" format library into an equivalent ENVI spectral library.  Extra header fields contain critical ancillary information like the record numbers and names. Syntax:
 
> specpr2envi [specpr_library_path] [output_envi_library_path] 

A header file will be created automatically. Then run the aggregation script.  It requires a translation file in comma-separated value format.  The translation file holds record numbers and spectrum names together with a 10-dimensional vector of mineral fraction coefficients for each spectrum.  These coefficients indicate how to portion the baand depth estimate into EMIT mineral classes.  An example of this format is: "convert_tetracorder/mineral_fractions_s06an14a.csv".  The aggregator syntax is: 

> python3 aggregator.py [tetra_expert_system_path] [tetra_enviFormat_library_path] [tetra_output_base_directory] [translation_file] [output_l2b]
