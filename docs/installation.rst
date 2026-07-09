Installation
============

Prerequisites
-------------

Before using RTTOVpy, the following software and data are required:

* RTTOV, preferably version 14 or later
* RTTOV coefficient files for the selected satellite and sensor
* Python 3
* WRF output files or ERA5 input data
* A Space-Track account if historical TLE retrieval is required

RTTOV must be downloaded and installed separately from the EUMETSAT NWP-SAF
website. RTTOVpy does not include RTTOV itself or RTTOV coefficient files.

Python Dependencies
-------------------

Core Python dependencies include:

.. list-table::
   :header-rows: 1

   * - Package
     - Purpose
   * - ``numpy``
     - Array operations
   * - ``netCDF4``
     - Reading WRF NetCDF files
   * - ``xarray``
     - NetCDF post-processing
   * - ``pyyaml``
     - Reading YAML namelists
   * - ``wrf-python``
     - WRF variable extraction
   * - ``pyorbital``
     - Satellite position and viewing geometry from TLE data
   * - ``requests``
     - TLE retrieval from CelesTrak or Space-Track
   * - ``matplotlib``
     - Plotting
   * - ``cartopy``
     - Geo-referenced map plots

Additional dependencies for verification are:

.. list-table::
   :header-rows: 1

   * - Package
     - Purpose
   * - ``satpy``
     - Reading satellite data files
   * - ``xesmf``
     - Regridding satellite data
   * - ``pyresample``
     - Area definitions for satellite swaths
   * - ``pyproj``
     - Coordinate reference system transformations
   * - ``geocat-viz``
     - Taylor diagrams

Clone the Repository
--------------------

Clone RTTOVpy from GitHub::

    git clone https://github.com/anikfal/rttovpy.git
    cd rttovpy
