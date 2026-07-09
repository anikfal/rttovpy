Satellite Data Readers
======================

The verification workflow reads satellite observations using Satpy. The
satellite reader is selected by ``satellite_sensor_id`` in the namelist.

Selected commonly used readers are shown below.

.. list-table::
   :header-rows: 1

   * - ID
     - Satellite / Sensor
     - Format
   * - 1
     - GOES-R ABI
     - Level 1b NetCDF
   * - 5
     - Himawari AHI
     - Level 1 HRIT
   * - 6
     - Himawari AHI
     - Level 1b HSD
   * - 11
     - S-NPP / JPSS ATMS
     - Level 1B NetCDF
   * - 12
     - NOAA 15--19 / Metop AVHRR
     - AAPP format
   * - 21
     - MTG FCI
     - Level-1c NetCDF
   * - 30
     - GOES Imager
     - Level 1 GOES-13/14/15 NetCDF
   * - 36
     - Sentinel-2 MSI
     - Level 1C SAFE
   * - 45
     - Sentinel-3 OLCI
     - Level 1B NetCDF
   * - 47
     - Landsat-8/9 OLI/TIRS
     - L1 GeoTIFF
   * - 52
     - MSG SEVIRI
     - Level 1b HRIT
   * - 53
     - MSG SEVIRI
     - Native format
   * - 59
     - JPSS VIIRS
     - Level 1b NetCDF
   * - 62
     - JPSS VIIRS
     - HDF5 SDR format

The complete list of supported readers is provided in
``wrf_data/satellite_sensors_id.yaml``.
