Workflow
========

RTTOVpy is organized as a staged workflow. The exact namelist options differ
between the ERA5 and WRF applications, but the general logic is the same.

Phase 1: Generate RTTOV Profiles
--------------------------------

In the first phase, RTTOVpy reads the meteorological input data and generates
RTTOV-compatible profile files.

For WRF, RTTOVpy extracts variables such as temperature, water vapour, 2 m
temperature, 10 m winds, skin temperature, terrain height, cloud fraction, and
surface masks. It also computes the half-level pressure grid required by
RTTOV.

For ERA5, RTTOVpy downloads or reads the required pressure-level and surface
variables, extracts the atmospheric and surface state, and creates one profile
for each ERA5 grid point.

In both workflows, RTTOVpy computes satellite and solar viewing geometry for
each grid point and writes ``prof-NNNNNN.dat`` files.

Phase 2: Run RTTOV
------------------

After the profile files are generated, RTTOVpy creates a ready-to-run shell
script for the RTTOV forward model.

Typical commands are::

    ./run_era5_example_fwd.sh ARCH=gfortran

or::

    ./run_wrf_example_fwd.sh ARCH=gfortran

The ``ARCH`` value must match the compiler architecture used when RTTOV was
compiled.

Phase 3: Post-processing
------------------------

After RTTOV has produced one text output file per input profile, RTTOVpy can
extract the simulated quantities and store them as NetCDF files.

Typical post-processed variables include:

* brightness temperature,
* radiance,
* overcast radiance,
* surface-to-space transmittance, and
* surface emissivity.

RTTOVpy can also generate PNG plots for each simulated band.

Phase 4: Verification
---------------------

For WRF simulations, RTTOVpy can optionally compare simulated radiances with
real satellite observations. Satellite data are read using Satpy, regridded to
the WRF grid, and compared statistically.

The verification step can generate Taylor diagrams and text-based statistics.
