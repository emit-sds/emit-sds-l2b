<h1 align="center"> emit-sds-l2b </h1>

_NOTE - See the **develop** branch - set as default - for the latest updates._

Welcome to the EMIT Level Level 2B science data system repository.  To understand how this repository is linked to the rest of the emit-sds repositories, please see [the repository guide](https://github.com/emit-sds/emit-main/wiki/Repository-Guide).

The repository supports two methods for creating mineral maps.  The first method runs and post-processes the output of tetracorder, aggregating the band depth maps from individual library spectra into a combined 10-channel product that shows the band depths of the 10 EMIT mineral classes. This occurs in two main steps: 

1) Generating mineral maps via tetracorder
2) Post-processing these maps to the EMIT-10 estimated spectral abundance

An example execution of these steps can be found in the run_tetracorder.sh script, though critically the tetracorder environment must be set up appropriatey first.  See the [Tetracorder](https://github.com/emit-sds/spectroscopy-tetracorder) repository for more info.  An example call using run_tetracorder.sh might look like:

> sh run_tetracorder.sh [tetracorder_output_directory] [aggregated_output_base] [l2a_reflectance_file] [l2a_reflectance_uncertainty_file] [emit_utils_path] [emit_sds_l2b_path] [tetracorder_cmd_base_path]

Ultimately, tetracorder is actually executed within the [emit-main](https://github.com/emit-sds/emit-main) repo directly without using run_tetracorder.sh.

Before running, compile the program "specpr2envi" (navigate to Spectral-Library-Reader-master and type "make") and use it to convert the tetracorder "specpr" format library into an equivalent ENVI spectral library.  Extra header fields contain critical ancillary information like the record numbers and names.  This step must be repeated for each reference library used (multiple are supported). Syntax:
> specpr2envi [specpr_library_path] [output_envi_library_path] 


