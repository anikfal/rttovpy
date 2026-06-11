# RTTOVpy

**A Python utility to run the [RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov/) radiative transfer model using WRF model outputs or ERA5 reanalysis data.**

RTTOVpy automates the full workflow: extracting atmospheric profiles from model output, formatting them for RTTOV, generating and executing the RTTOV driver script, postprocessing the outputs into NetCDF files and publication-quality maps, and optionally verifying simulated radiances against real satellite observations.

---

## Table of Contents

1. [What is RTTOVpy?](#what-is-rttovpy)
2. [Features](#features)
3. [Prerequisites](#prerequisites)
4. [Python Dependencies](#python-dependencies)
5. [Repository Structure](#repository-structure)
6. [Configuration — `namelist_wrf.yaml`](#configuration--namelist_wrfyaml)
7. [Workflow Overview](#workflow-overview)
8. [Usage](#usage)
   - [Inspect a coefficient file](#inspect-a-coefficient-file)
   - [Inspect a WRF file](#inspect-a-wrf-file)
   - [Phase 1 — Generate profiles and run RTTOV](#phase-1--generate-profiles-and-run-rttov)
   - [Phase 2 — Postprocessing](#phase-2--postprocessing)
   - [Phase 3 — Verification](#phase-3--verification)
9. [Satellite Platforms](#satellite-platforms)
10. [Satellite Data Readers (Verification)](#satellite-data-readers-verification)
11. [Output Description](#output-description)
12. [Examples](#examples)
13. [Notes and Limitations](#notes-and-limitations)
14. [Author](#author)

---

## What is RTTOVpy?

[RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov/) (Radiative Transfer for TOVS) is a fast radiative transfer model developed by EUMETSAT NWP-SAF, widely used for satellite data assimilation and forward simulation of satellite observations from atmospheric model output. RTTOV itself is a Fortran library that requires per-profile ASCII input files and a Fortran/C driver program.

RTTOVpy bridges the gap between Python-based atmospheric model data (WRF NetCDF, ERA5 GRIB) and RTTOV's Fortran world. Given a WRF output file and a short YAML configuration, RTTOVpy will:

- Extract every atmospheric column in the WRF domain as a properly formatted RTTOV profile
- Compute accurate satellite viewing geometry from orbital elements (TLE)
- Generate and run the RTTOV driver shell script
- Collect the outputs into geo-referenced NetCDF files and maps
- Optionally verify the simulated radiances against real satellite observations

---

## Features

- **WRF and ERA5 support** — two parallel pipelines, each with its own namelist
- **RTTOV 14+ and legacy 12/13 support** — automatically dispatched from a single entry point
- **Automatic satellite geometry** — zenith and azimuth angles computed from TLE data; no manual angle input required
- **Dual TLE sources** — uses [CelesTrak](https://celestrak.org/) for recent observations, [Space-Track](https://www.space-track.org/) for historical TLE retrieval (configurable)
- **Solar channel auto-detection** — inspects the RTTOV coefficient file to determine whether shortwave channels are selected, and enables solar simulation automatically
- **WRF-Chem dust / aerosol support** — ingests `DUST_1`–`DUST_5` from WRF-Chem output and generates RTTOV aerosol profile files
- **Postprocessing to NetCDF + maps** — produces brightness temperature, radiance, overcast radiances, surface-to-space transmittance, and emissivities as geo-referenced NetCDF files with Cartopy plots
- **RGB composite plots** — user-defined band arithmetic for RGB brightness temperature composites
- **Satellite data verification** — regrids real satellite data (via [satpy](https://satpy.readthedocs.io/) + [xesmf](https://xesmf.readthedocs.io/)) onto the WRF grid and generates Taylor diagrams
- **WRF projection-aware plotting** — Lambert Conformal, Mercator, and Polar Stereographic projections are detected automatically from WRF metadata

---

## Prerequisites

- **RTTOV** (version 14 or later recommended) — installed and compiled. Obtain from [NWP-SAF](https://nwp-saf.eumetsat.int/site/software/rttov/download/).
- **RTTOV coefficient files** — download for your specific satellite/sensor from the NWP-SAF website.
- **WRF output file** (NetCDF) — or ERA5 GRIB files for the ERA5 pipeline.
- A [Space-Track](https://www.space-track.org/) account if you need historical TLE data (observations older than ~2 days).

---

## Python Dependencies

Install via `pip` or `conda`. Core dependencies:

| Package | Purpose |
|---|---|
| `numpy` | Array operations |
| `netCDF4` | Reading WRF files |
| `xarray` | NetCDF postprocessing |
| `pyyaml` | Namelist parsing |
| `wrf-python` | WRF variable extraction |
| `pyorbital` | Satellite orbital calculations from TLE |
| `requests` | TLE retrieval from CelesTrak / Space-Track |
| `matplotlib` | Plotting |
| `cartopy` | Geo-referenced map plots |

Additional dependencies for verification:

| Package | Purpose |
|---|---|
| `satpy` | Reading satellite data files (62+ formats) |
| `xesmf` | Regridding satellite swath to WRF grid |
| `pyresample` | Area definitions for regridding |
| `pyproj` | CRS transformations |
| `geocat-viz` | Taylor diagram |

---

## Repository Structure

```
rttovpy/
├── README.md
│
├── wrf_data/                          # WRF pipeline
│   ├── rttovpy.py                     # Entry point
│   ├── namelist_wrf.yaml              # User configuration
│   ├── satellite_names.yaml           # Satellite platform index
│   ├── satellite_sensors_id.yaml      # Satpy reader index (for verification)
│   │
│   ├── modules14plus/                 # RTTOV 14+ modules
│   │   ├── run_rttov14.py             # Main pipeline logic
│   │   ├── p_stag.py                  # WRF staggered pressure computation
│   │   ├── rttov_utils.py             # Coefficient file inspection, domain/time helpers
│   │   ├── application_shell.py       # RTTOV driver shell script generator
│   │   ├── tle_fetcher.py             # TLE retrieval (CelesTrak + Space-Track)
│   │   ├── choose_tle_source.py       # Selects TLE source based on observation age
│   │   ├── satpy_readers.py           # Satpy reader ID → reader name mapping
│   │   ├── satellite_celestrak_urls.yaml   # CelesTrak URLs per satellite
│   │   ├── satellite_to_norad_id.yaml      # NORAD IDs for Space-Track queries
│   │   ├── run_wrf_example_fwd.sh     # Shell script template (standard)
│   │   └── run_wrfchem_dust_example_fwd.sh # Shell script template (dust/aerosol)
│   │
│   └── modules/                       # Legacy modules for RTTOV 12/13
│
└── era5_input_data/                   # ERA5 pipeline
    ├── rttovpy.py
    └── namelist_era5.yaml
```

---

## Configuration — `namelist_wrf.yaml`

All settings are controlled through `wrf_data/namelist_wrf.yaml`. The key sections are:

### RTTOV paths

```yaml
rttov_version: 14
rttov_installation_path: /path/to/rttov14
rttov_coefficient_file_path: /path/to/rttov14/rtcoef_rttov14/.../rtcoef_<platform>_<sensor>.dat
```

The coefficient file must match the satellite platform selected below. RTTOVpy validates this automatically and will exit with a clear error if they don't match.

### Input data

```yaml
wrf_file_path: /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00
```

### Simulation time

```yaml
time_of_simulation:
  year: 2022
  month: 3
  day: 13
  hour: 12        # UTC hour; must fall within the time range of the WRF file
```

### Solar simulation

```yaml
solar_simulation:
  enabled: true   # Required for shortwave (< 4 µm) channels
                  # Auto-enabled if solar channels are detected in the coefficient file
```

### Satellite settings

```yaml
satellite_information:
  sat_name_index: 20         # Index from satellite_names.yaml (1–53; see full list below)
  sat_channel_list: [17, 18, 19]         # RTTOV channel numbers in the coefficient file
  sat_channel_names: [M12, M13, M14]    # Corresponding band names in the satellite data file

  user_defined_position:
    enabled: true            # Use a fixed satellite position (overrides TLE-computed position)
    sat_latitude: 33
    sat_longitude: 51

  historical_tle:
    enabled: true            # Use Space-Track for observations older than ~2 days
    space-track.org_username: your_username
    space-track.org_password: your_password
```

If `user_defined_position` is disabled, RTTOVpy computes the satellite position from TLE data at the exact observation time.

### WRF-Chem dust / aerosol

```yaml
wrfchem_dust_profiles:
  enabled: false
  aerosol_coefficient_file_path: /path/to/rttov14/rtcoef_rttov14/aertable_visir/rttov_aertable_<platform>_<sensor>.dat
```

When enabled, RTTOVpy reads `DUST_1`–`DUST_5` variables from WRF-Chem output and maps them to RTTOV aerosol types (mineral nucleation, accumulation, coarse modes).

### Postprocessing

```yaml
postprocessing:
  enabled: true
  postprocessing_directory_suffix: rttov_outputs_postprocessing
  image_plot_all_bands: true
  RGB_plot_brightness_temperature:
    enabled: false
    Red: B5 + B6 + B7    # Band arithmetic expressions
    Green: B5
    Blue: B7
```

### Verification

```yaml
verification:
  enabled: true
  verification_directory_suffix: viirs_verification
  satellite_file_path: /path/to/satellite_file.h5
  satellite_files_group:
    enabled: true              # Load all files in a directory instead of a single file
    satellite_file_directory: /path/to/satellite_files/
  satellite_sensor_id: 62     # Index from satellite_sensors_id.yaml (1–62; see full list below)
  taylor_diagram_name: radiation_taylor_diagram
  keep_remapped_satellite_to_wrf_data:
    enabled: true
    remapped_file_name: viirs_to_wrf   # Prefix for the output remapped NetCDF file
```

---

## Workflow Overview

The pipeline has three phases, each controlled by a flag in `namelist_wrf.yaml`. Typically you run each phase sequentially in separate calls:

```
Phase 1: make_inputdata()
    WRF NetCDF → per-column prof-NNNNNN.dat files → run_wrf_example_fwd.sh

    [Run RTTOV manually: ./run_wrf_example_fwd.sh ARCH=gfortran]

Phase 2: run_postprocessing()  (postprocessing.enabled = true)
    RTTOV text outputs → NetCDF files + PNG maps

Phase 3: verification()  (verification.enabled = true)
    Real satellite data + RTTOV radiances → Taylor diagram + statistics
```

### Phase 1 — Profile extraction and RTTOV script generation

For each WRF grid column, RTTOVpy:

1. Reads the WRF variables: temperature (`tk`), water vapour (`QVAPOR`), 2m temperature (`T2`), 10m winds (`U10`, `V10`), skin temperature (`TSK`), terrain height (`HGT`), cloud fraction (`cloudfrac`), lake mask (`LAKEMASK`).
2. Computes staggered (half-level) pressures ordered TOA → surface using height interpolation (`p_stag.py`). This maps WRF's full-level pressure grid to the half-level grid required by RTTOV.
3. Fetches the satellite TLE and computes the satellite's position at the observation time (`tle_fetcher.py`).
4. Computes per-column satellite zenith/azimuth and solar zenith/azimuth angles using `pyorbital`.
5. Writes a `prof-NNNNNN.dat` file for each column in RTTOV's ASCII profile format.
6. Generates `run_wrf_example_fwd.sh` — a ready-to-run shell script that calls RTTOV's `example_fwd.exe` for every profile.

After Phase 1, run RTTOV:

```bash
./run_wrf_example_fwd.sh ARCH=gfortran
```

This produces one output file per profile in the `*_outputs222/` directory.

### Phase 2 — Postprocessing

Reads all RTTOV output files, extracts:

- Brightness temperatures (K)
- Calculated radiances (mW/m²/sr/cm⁻¹)
- Overcast radiances
- Surface-to-space transmittance
- Surface emissivities

Saves each variable as a geo-referenced xarray/NetCDF Dataset retaining the WRF projection metadata (Lambert Conformal, Mercator, Polar Stereographic). Generates Cartopy maps for each band.

### Phase 3 — Verification

Loads real satellite data using satpy (see [Satellite Data Readers](#satellite-data-readers-verification)), regrids the swath onto the WRF domain using bilinear interpolation with xesmf, then computes:

| Statistic | Description |
|---|---|
| CV (coefficient of variation of standard deviation) | Normalised spread relative to observations |
| RMSE | Normalised root-mean-square error |
| Correlation | Pearson correlation coefficient |

Outputs a Taylor diagram (PNG) and a statistics table (`.txt`) to the verification directory.

---

## Usage

All commands are run from `wrf_data/`:

### Inspect a coefficient file

List all channels with their wavenumber, wavelength, and whether solar simulation is required:

```bash
python rttovpy.py --wavelength /path/to/rtcoef_<platform>_<sensor>.dat
```

Example output:
```
Channel | Wavenumber (cm⁻¹) | Wavelength (µm) | Solar
------------------------------------------------------------
      1 |          748.000  |          13.369 | False
     17 |         2662.000  |           3.757 | False
     18 |         1000.000  |          10.000 | False
```

### Inspect a WRF file

Print the time range of a WRF output file:

```bash
python rttovpy.py --wrftime /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00
```

Print the lat/lon domain extent and midpoint:

```bash
python rttovpy.py --wrfdomain /path/to/wrfout_d01_YYYY-MM-DD_HH:00:00
```

### Phase 1 — Generate profiles and run RTTOV

Set `postprocessing.enabled: false` and `verification.enabled: false` in the namelist, then:

```bash
python rttovpy.py
./run_wrf_example_fwd.sh ARCH=gfortran
```

If the profile directory already exists (from a previous run), RTTOVpy will skip re-extraction and only regenerate the shell script — useful when you want to change channels or coefficient files without reprocessing all profiles.

### Phase 2 — Postprocessing

Set `postprocessing.enabled: true` and `verification.enabled: false`:

```bash
python rttovpy.py
```

### Phase 3 — Verification

Set `verification.enabled: true`. This phase also requires Phase 2 outputs (the radiance NetCDF) to already exist:

```bash
python rttovpy.py
```

---

## Satellite Platforms

The satellite platform is selected by `sat_name_index` in the namelist. The full list (`satellite_names.yaml`):

| Index | Platform | Index | Platform | Index | Platform |
|---|---|---|---|---|---|
| 1 | EOS-Aqua | 19 | NOAA-19 | 37 | Metop-B |
| 2 | EOS-Terra | 20 | NOAA-20 | 38 | Metop-C |
| 3 | Landsat-7 | 21 | Meteosat-7 | 39 | SMOS |
| 4 | Landsat-8 | 22 | Meteosat-8 | 40 | CloudSat |
| 5 | Sentinel-3A | 23 | Meteosat-9 | 41 | CALIPSO |
| 6 | Sentinel-3B | 24 | Meteosat-10 | 42 | DMSP-F15 |
| 7 | NOAA-6 | 25 | Meteosat-11 | 43 | DMSP-F16 |
| 8 | NOAA-7 | 26 | Meteosat-12 | 44 | DMSP-F17 |
| 9 | NOAA-8 | 27 | GOES-13 | 45 | DMSP-F18 |
| 10 | NOAA-9 | 28 | GOES-14 | 46 | DMSP-F19 |
| 11 | NOAA-10 | 29 | GOES-15 | 47 | FY-2D |
| 12 | NOAA-11 | 30 | GOES-16 | 48 | FY-2E |
| 13 | NOAA-12 | 31 | Suomi-NPP | 49 | FY-2G |
| 14 | NOAA-14 | 32 | Himawari-8 | 50 | FY-3A |
| 15 | NOAA-15 | 33 | Himawari-9 | 51 | FY-3B |
| 16 | NOAA-16 | 34 | INSAT-3D | 52 | FY-3C |
| 17 | NOAA-17 | 35 | JASON-2 | 53 | FY-3D |
| 18 | NOAA-18 | 36 | Metop-A | | |

The satellite name must be consistent with the RTTOV coefficient file name. RTTOVpy validates this at startup.

---

## Satellite Data Readers (Verification)

The verification sensor is selected by `satellite_sensor_id` in the namelist. These map to satpy reader names. Selected commonly-used readers:

| ID | Satellite / Sensor | Format |
|---|---|---|
| 1 | GOES-R ABI | Level 1b NetCDF |
| 5 | Himawari AHI | Level 1 HRIT |
| 6 | Himawari AHI | Level 1b HSD |
| 11 | S-NPP / JPSS ATMS | Level 1B NetCDF |
| 12 | NOAA 15–19 / Metop AVHRR | AAPP format |
| 21 | MTG FCI | Level-1c NetCDF |
| 30 | GOES Imager | Level 1 GOES-13/14/15 netCDF |
| 36 | Sentinel-2 MSI | Level 1C SAFE |
| 45 | Sentinel-3 OLCI | Level 1B NetCDF |
| 47 | Landsat-8/9 OLI/TIRS | L1 GeoTIFF |
| 52 | MSG SEVIRI | Level 1b HRIT |
| 53 | MSG SEVIRI | Native format |
| 59 | JPSS VIIRS | Level 1b NetCDF |
| 62 | JPSS VIIRS | HDF5 SDR format |

The full list of all 62 supported readers is in `wrf_data/satellite_sensors_id.yaml`.

---

## Output Description

RTTOVpy creates the following directories alongside the WRF output file:

| Directory | Contents |
|---|---|
| `<wrffile>_inputs222/` | Per-column RTTOV profile files (`prof-NNNNNN.dat`) |
| `<wrffile>_outputs222/` | Per-column RTTOV output files (created by RTTOV, read by postprocessing) |
| `<wrffile>_rttov_outputs_postprocessing/` | NetCDF files and PNG maps for each output variable and band |
| `<wrffile>_<verification_suffix>/` | Remapped satellite NetCDF files, Taylor diagram PNG, statistics table |

NetCDF files preserve WRF projection metadata (central lat/lon, standard parallels, projection type) so they can be opened in any GIS tool.

---

## Examples

> *Examples will be added here.*

---

## Notes and Limitations

- **Profile count**: RTTOVpy processes every grid column in the WRF domain, which can be time-consuming for large domains. RTTOV is called once per profile via the shell script.
- **Zenith angle cutoff**: RTTOV skips profiles where the satellite zenith angle exceeds ~75°. Those output files are filled with `missing_value`.
- **Pressure levels**: RTTOVpy computes RTTOV's required half-level pressure grid from WRF's staggered geopotential heights using log-linear interpolation. The surface level is taken from PSFC directly, and the top level is extrapolated.
- **WRF-Chem aerosol**: RTTOV maps WRF-Chem dust bins to three aerosol types: mineral nucleation (DUST_1), mineral accumulation (DUST_2 + DUST_3 + DUST_4), and mineral coarse (DUST_5).
- **TLE accuracy**: For historical simulations (observations older than ~2 days), Space-Track credentials are required. For recent observations, CelesTrak is used automatically.
- **Coefficient file / satellite consistency**: The satellite name (from `sat_name_index`) must match the coefficient file path. RTTOVpy checks this and exits with a descriptive error if they don't match.

---

## Author

**Amirhossein Nikfal**
[GitHub](https://github.com/anikfal)
