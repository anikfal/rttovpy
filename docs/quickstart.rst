Quick Start
===========

ERA5 Workflow
-------------

A basic ERA5-based simulation can be run from the ``era5_input_data`` directory::

    cd era5_input_data
    python rttovpy.py
    ./run_era5_example_fwd.sh ARCH=gfortran

After the RTTOV forward model finishes, enable post-processing in
``namelist_era5.yaml`` and run RTTOVpy again::

    python rttovpy.py

WRF Workflow
------------

A basic WRF-based simulation can be run from the ``wrf_data`` directory::

    cd wrf_data
    python rttovpy.py
    ./run_wrf_example_fwd.sh ARCH=gfortran

After the RTTOV forward model finishes, enable post-processing in
``namelist_wrf.yaml`` and run RTTOVpy again::

    python rttovpy.py

Inspect a Coefficient File
--------------------------

To list the channels available in an RTTOV coefficient file, run::

    python rttovpy.py --wavelength /path/to/rtcoef_<platform>_<sensor>.dat

This prints the RTTOV channel number, central wavenumber, wavelength, and
whether the channel requires solar simulation.

Inspect a WRF File
------------------

To print the time range of a WRF output file::

    python rttovpy.py --wrftime /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00

To print the latitude/longitude domain extent and midpoint::

    python rttovpy.py --wrfdomain /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00
