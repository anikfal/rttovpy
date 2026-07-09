WRF Configuration
=================

The WRF workflow is controlled by ``wrf_data/namelist_wrf.yaml``.

RTTOV Paths
-----------

.. code-block:: yaml

    rttov_version: 14
    rttov_installation_path: /path/to/rttov14
    rttov_coefficient_file_path: /path/to/rtcoef_<platform>_<sensor>.dat

``rttov_version``
    RTTOV major version.

``rttov_installation_path``
    Root directory of the RTTOV installation.

``rttov_coefficient_file_path``
    RTTOV coefficient file corresponding to the selected satellite and sensor.

Input Data
----------

.. code-block:: yaml

    wrf_file_path: /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00

This option defines the WRF output file used to generate RTTOV profiles.

Simulation Time
---------------

.. code-block:: yaml

    time_of_simulation:
      year: 2022
      month: 3
      day: 13
      hour: 12

The selected time must fall within the time range available in the WRF output
file.

Solar Simulation
----------------

.. code-block:: yaml

    solar_simulation:
      enabled: true

Solar simulation is required for shortwave channels. RTTOVpy can inspect the
selected coefficient file and automatically enable solar simulation when solar
channels are detected.

Satellite Settings
------------------

.. code-block:: yaml

    satellite_information:
      sat_name_index: 20
      sat_channel_list: [17, 18, 19]
      sat_channel_names: [M12, M13, M14]

      user_defined_position:
        enabled: true
        sat_latitude: 33
        sat_longitude: 51

      historical_tle:
        enabled: true
        space-track.org_username: your_username
        space-track.org_password: your_password

``sat_name_index``
    Index of the selected satellite platform from ``satellite_names.yaml``.

``sat_channel_list``
    RTTOV channel numbers to simulate.

``sat_channel_names``
    User-defined channel labels used in the output files.

``user_defined_position``
    Uses a fixed satellite position instead of a TLE-derived position.

``historical_tle``
    Enables historical TLE retrieval from Space-Track.

WRF-Chem Dust and Aerosol
-------------------------

.. code-block:: yaml

    wrfchem_dust_profiles:
      enabled: false
      aerosol_coefficient_file_path: /path/to/rttov_aertable_<platform>_<sensor>.dat

When enabled, RTTOVpy reads WRF-Chem dust variables and maps them to RTTOV
aerosol profile inputs.

Post-processing
---------------

.. code-block:: yaml

    postprocessing:
      enabled: true
      postprocessing_directory_suffix: rttov_outputs_postprocessing
      image_plot_all_bands: true
      RGB_plot_brightness_temperature:
        enabled: false
        Red: B5 + B6 + B7
        Green: B5
        Blue: B7

Verification
------------

.. code-block:: yaml

    verification:
      enabled: true
      verification_directory_suffix: viirs_verification
      satellite_file_path: /path/to/satellite_file.h5
      satellite_files_group:
        enabled: true
        satellite_file_directory: /path/to/satellite_files/
      satellite_sensor_id: 62
      taylor_diagram_name: radiation_taylor_diagram
      keep_remapped_satellite_to_wrf_data:
        enabled: true
        remapped_file_name: viirs_to_wrf
