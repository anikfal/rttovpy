Overview
========

What is RTTOVpy?
----------------

RTTOVpy is a Python utility for running the RTTOV radiative transfer model
using WRF model outputs or ERA5 reanalysis data.

RTTOV itself is a Fortran radiative transfer model that requires atmospheric
profiles, surface information, satellite viewing geometry, coefficient files,
and a forward-model driver program. RTTOVpy connects Python-based atmospheric
data workflows with RTTOV by preparing the required input profiles and
organizing the model outputs.

Main Capabilities
-----------------

RTTOVpy can

* extract atmospheric and surface profiles from WRF output,
* download and process ERA5 reanalysis data,
* generate RTTOV-compatible ASCII profile files,
* compute satellite viewing geometry from TLE data,
* generate shell scripts for RTTOV forward simulations,
* convert RTTOV text outputs into NetCDF files,
* generate quick-look satellite image products,
* create RGB composites from user-defined band expressions, and
* optionally verify simulated radiances against real satellite observations.

Supported Workflows
-------------------

RTTOVpy currently provides two main workflows:

``era5_input_data/``
    ERA5-based simulations. This workflow downloads ERA5 data, extracts the
    required RTTOV input variables, runs RTTOV, and post-processes the
    simulated satellite quantities.

``wrf_data/``
    WRF-based simulations. This workflow extracts atmospheric profiles from
    WRF or WRF-Chem output and can also perform verification against real
    satellite observations.

Repository Layout
-----------------

The main repository structure is::

    rttovpy/
    ├── docs/
    ├── era5_input_data/
    ├── wrf_data/
    └── README.md

The ``README.md`` file provides a short project introduction, while this
documentation contains the detailed user guide and examples.
