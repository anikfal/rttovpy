Notes and Limitations
=====================

Profile Count
-------------

RTTOVpy processes every grid column in the selected WRF or ERA5 domain. This
can be time-consuming for large domains because RTTOV is called once per
profile by the generated shell script.

Zenith Angle Cutoff
-------------------

RTTOVpy skips profiles where the satellite zenith angle exceeds the supported
limit. These locations are filled with missing values in the post-processed
outputs.

Pressure Levels
---------------

For WRF input, RTTOVpy computes the required half-level pressure grid from WRF
model data. The pressure grid is ordered from the top of the atmosphere to the
surface, following the RTTOV profile format.

WRF-Chem Aerosol
----------------

For WRF-Chem simulations, RTTOVpy maps WRF-Chem dust bins to RTTOV aerosol
types. The current dust mapping is:

* ``DUST_1`` to mineral nucleation mode,
* ``DUST_2`` + ``DUST_3`` + ``DUST_4`` to mineral accumulation mode,
* ``DUST_5`` to mineral coarse mode.

TLE Accuracy
------------

For recent simulations, RTTOVpy can use CelesTrak TLE data. For historical
simulations or higher-precision orbital information, Space-Track credentials
can be provided in the namelist.

Coefficient File Consistency
----------------------------

The satellite selected by ``sat_name_index`` must be consistent with the RTTOV
coefficient file path. RTTOVpy checks this consistency and reports a clear
message if the selected satellite and coefficient file do not match.
