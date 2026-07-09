Satellite Platforms
===================

The satellite platform is selected by ``sat_name_index`` in the namelist. The
index is defined in ``satellite_names.yaml``.

.. list-table::
   :header-rows: 1

   * - Index
     - Platform
     - Index
     - Platform
     - Index
     - Platform
   * - 1
     - EOS-Aqua
     - 19
     - NOAA-19
     - 37
     - Metop-B
   * - 2
     - EOS-Terra
     - 20
     - NOAA-20
     - 38
     - Metop-C
   * - 3
     - Landsat-7
     - 21
     - Meteosat-7
     - 39
     - SMOS
   * - 4
     - Landsat-8
     - 22
     - Meteosat-8
     - 40
     - CloudSat
   * - 5
     - Sentinel-3A
     - 23
     - Meteosat-9
     - 41
     - CALIPSO
   * - 6
     - Sentinel-3B
     - 24
     - Meteosat-10
     - 42
     - DMSP-F15
   * - 7
     - NOAA-6
     - 25
     - Meteosat-11
     - 43
     - DMSP-F16
   * - 8
     - NOAA-7
     - 26
     - Meteosat-12
     - 44
     - DMSP-F17
   * - 9
     - NOAA-8
     - 27
     - GOES-13
     - 45
     - DMSP-F18
   * - 10
     - NOAA-9
     - 28
     - GOES-14
     - 46
     - DMSP-F19
   * - 11
     - NOAA-10
     - 29
     - GOES-15
     - 47
     - FY-2D
   * - 12
     - NOAA-11
     - 30
     - GOES-16
     - 48
     - FY-2E
   * - 13
     - NOAA-12
     - 31
     - Suomi-NPP
     - 49
     - FY-2G
   * - 14
     - NOAA-14
     - 32
     - Himawari-8
     - 50
     - FY-3A
   * - 15
     - NOAA-15
     - 33
     - Himawari-9
     - 51
     - FY-3B
   * - 16
     - NOAA-16
     - 34
     - INSAT-3D
     - 52
     - FY-3C
   * - 17
     - NOAA-17
     - 35
     - JASON-2
     - 53
     - FY-3D
   * - 18
     - NOAA-18
     - 36
     - Metop-A
     -
     -

The satellite name must be consistent with the selected RTTOV coefficient file.
RTTOVpy checks this consistency at startup.
