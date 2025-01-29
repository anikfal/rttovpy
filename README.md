# RTTOVpy

**A python utility to run the [RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov/) model, using WRF model outputs or ERA5 data.**

Easily run the RTTOV model by setting up a namelist file.

## Running the RTTOV model with ERA5 data:

-    Specify the simulation time.
-    Specify the geographical domain of the simulation.
-    Specify the satellite sensor and its specifications.
-    Download the required ERA5 data for running the RTTOV model.
-    Create a shell script to run the RTTOV model using the downloaded ERA5 data.
-    Run the RTTOV model using the generated shell script for each grid point in the downloaded ERA5 data.
-    Create a NetCDF file from the RTTOV outputs for each grid point.

## Running the RTTOV model with WRF model output:

-    Specify the simulation time.
-    Specify the satellite sensor and its specifications.
-    Extract the required atmospheric data from the WRF model output for running the RTTOV model.
-    Create a shell script to run the RTTOV model using the extracted data from the WRF output.
-    Run the RTTOV model using the generated shell script for each grid point in the WRF model outputs.
-    Create a NetCDF file from the RTTOV outputs for each grid point.
-    Plot the NetCDF data for each band.
-    Create RGB images of brightness temperature for selected bands.