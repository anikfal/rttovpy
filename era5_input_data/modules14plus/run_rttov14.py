## Code for downloading ERA5 data,and making RTTOV profile data based on the downloaded data,
## As well as generating the shellscript for running the RTTOV model using the profile data
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################
import importlib
required_modules = ["numpy", "netCDF4", "pyorbital", "cfgrib", "xarray", "yaml"]
for module in required_modules:
    try:
        importlib.import_module(module)
    except:
        print("Warning: The Python module", module, "is not installed.")
        print("Install it and run again.")
        print("Exiting ..")
        exit()

from modules14plus import era5_download_manager, count_lines, application_shell
from modules14plus.conversions import surface_humidity
from modules14plus.rttov_utils import requires_solar_for_channels
import yaml, os
from glob import glob
from pyorbital.orbital import get_observer_look
from pyorbital.orbital import Orbital
from pyorbital.astronomy import get_alt_az
from datetime import datetime, timezone
from math import pi

with open('namelist_era5.yaml', 'r') as yaml_file:
    namelist = yaml.safe_load(yaml_file)
satIndex = namelist["satellite_information"]["sat_name_index"]
satChannels000 = namelist["satellite_information"]["sat_channel_list"]
satChannels = " ".join(str(item_value) for item_value in satChannels000)
bandNames = namelist["satellite_information"]["sat_channel_names"]
if len(satChannels000) != len(bandNames):
    print("Warning: lengths of sat_channel_list and sat_channel_names are not equal.")
    print("Exiting ..")
    exit()
rttov_install_path = namelist["rttov_installation_path"]
if rttov_install_path.endswith('/'):
    rttov_install_path = rttov_install_path[:-1]
if satIndex > 76:
    print("Warning: Satellite name index must be between 1 to 76. Please look inside <satellite_names.yaml>.")
    print("exiting ..")
    exit()
angleEnable = namelist["satellite_information"]["user_defined_sat_position"]['enabled']

with open('satellite_names.yaml', 'r') as yaml_file:
    satNameFile = yaml.safe_load(yaml_file)

with open('modules14plus/satellite_celestrak_urls.yaml', 'r') as yaml_file:
    satCelestrakUrls = yaml.safe_load(yaml_file)

filePrefix = namelist["area_of_simulation"]["domain_name"]
dirName = filePrefix+"_profiles/"
era5_surface_file000 = glob("era5data_surface_level_*" + filePrefix + ".grib")#[0]
era5_level_file000 = glob("era5data_pressure_levels_*" + filePrefix + ".grib")#[0]
rttovCoef = namelist["rttov_coefficient_file_path"]
postprocessingEnabled = namelist["postprocessing"]['enabled']

import numpy as np
import xarray as xr

def make_inputdata():
    do_solar = int(namelist["solar_simulation"]['enabled'])
    # --- CHECK SOLAR REQUIREMENT ---
    need_solar, solar_channels, channel_info = requires_solar_for_channels(
        rttovCoef,
        satChannels000
    )
    if need_solar:
        print(f"Channels requiring solar: {solar_channels}")
        if not do_solar:
            print("WARNING: solar_simulation is disabled, but solar channels detected.")
            do_solar = 1
            print(" solar_simulation is automatically set to true.")
            # Option 2 (strict mode)
            # raise RuntimeError("Solar channels selected but DO_SOLAR=0")

    if os.path.exists(dirName) and os.path.isdir(dirName):
        print("Using existing profile directory:", dirName)
        profiles = [f for f in os.listdir(dirName) if f.startswith("prof-")]
        if profiles:
            pressureLevelsSize = count_lines.count_lines_between(
                dirName + profiles[0],
                "! Pressure levels (hPa)",
                "! Temperature profile (K)"
            ) - 2
            application_shell.make_final_application_shell(
                rttovCoef,
                str(pressureLevelsSize),
                str(do_solar),
                satChannels,
                str(len(satChannels000)),
                rttov_install_path
            )
            return

    print("Making profile files for RTTOV ..")
    if ( # check if ERA5 files are available before downloading
        len(era5_surface_file000) > 0 and
        len(era5_level_file000) > 0 and
        os.path.exists(era5_surface_file000[0]) and
        os.path.exists(era5_level_file000[0])
    ):
        era5_surface_file = era5_surface_file000[0]
        era5_level_file = era5_level_file000[0]

        print(f"Using existing ERA5 files:")
        print(f"  {era5_surface_file}")
        print(f"  {era5_level_file}")

    else:
        print("ERA5 files not found. Downloading...")
        era5_download_manager.main_dm()

        era5_surface_file = glob("era5data_surface_level_*" + filePrefix + ".grib")[0]
        era5_level_file = glob("era5data_pressure_levels_*" + filePrefix + ".grib")[0]
    
    os.makedirs(dirName)  # no need to check if it exists, as it's already been checked before

    level_ds = xr.open_dataset(
        era5_level_file,
        engine="cfgrib",
        backend_kwargs={"read_keys": ["pv", "NV"], "indexpath": ""}
    )

    surface_ds = xr.open_dataset(
        era5_surface_file,
        engine="cfgrib",
        backend_kwargs={"indexpath": ""}
    )

    temperature = level_ds.t.values
    qv = level_ds.q.values

    t2m = surface_ds.t2m.values
    d2m = surface_ds.d2m.values
    sp = surface_ds.sp.values
    u10 = surface_ds.u10.values
    v10 = surface_ds.v10.values
    skinT = surface_ds.skt.values
    landSeaMask = surface_ds.lsm.values
    geopotential = surface_ds.z.values
    cloudFraction = surface_ds.tcc.values

    lat = surface_ds.latitude.values
    lon = surface_ds.longitude.values

    # --- pressure ---
    pv = np.array(level_ds.t.attrs["GRIB_pv"])
    n_half = len(pv) // 2
    a_half = pv[:n_half]
    b_half = pv[n_half:]

    a_3d = a_half[:, None, None]
    b_3d = b_half[:, None, None]

    p_half = a_3d + b_3d * sp

    # --- FIX top level using log extrapolation ---
    p_half[0, :, :] = (p_half[1, :, :] ** 2) / p_half[2, :, :]

    q2m = surface_humidity(d2m, sp)

    year = namelist["time_of_simulation"]["year"]
    month = namelist["time_of_simulation"]["month"]
    day = namelist["time_of_simulation"]["day"]
    hour = namelist["time_of_simulation"]["hour"]

    observationTime = datetime(year, month, day, hour)


    # import requests
    # session = requests.Session()
    # session.headers.update({"User-Agent": "Mozilla/5.0"})
    # payload = {
    #     "identity": "ah.nikfal@gmail.com",
    #     "password": "kjdfkjdf1234PP!."  # remember to change this!
    # }
    # login = session.post("https://www.space-track.org/ajaxauth/login", data=payload)
    # print(login.status_code, login.text)



    import requests
    import tempfile
    celestrak_url = satCelestrakUrls[satNameFile[satIndex]]
    tel_text = requests.get(celestrak_url, timeout=10).text
    # write to temp file
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp.write(tel_text)
    tmp.close()
    try:
        orb = Orbital(satNameFile[satIndex], tle_file=tmp.name)
        print("Retrieving satellite specific information for", satNameFile[satIndex], "for the observation time:", observationTime)
        satPositions= orb.get_lonlatalt(observationTime) #Get longitude, latitude and altitude of the satellite
    finally:
        os.remove(tmp.name)



    # orb = Orbital(satNameFile[satIndex])
    # orb = Orbital(satNameFile[satIndex], tle_file="...")
    # print("Retrieving satellite specific information for", satNameFile[satIndex], "for the observation time:", observationTime)
    # print("please wait ..")
    # satPositions= orb.get_lonlatalt(observationTime) #Get longitude, latitude and altitude of the satellite
    # satPositions = (136.85902196460546, -53.70781534686423, 715.6113205704698)
    satAltitude = satPositions[2]

    if angleEnable:
        satLat = namelist["satellite_information"]["user_defined_sat_position"]["sat_latitude"]
        satLon = namelist["satellite_information"]["user_defined_sat_position"]["sat_longitude"]
    else:
        satLat = satPositions[1]
        satLon = satPositions[0]

    nlev, jjmax, iimax = temperature.shape
    # print("Number of levels:", nlev+1)
    # exit()
    profileCount = 1
    header_text = """! Specify input profiles for example_fwd.F90 and other example programs.
    ! Multiple profiles may be described: follow the same format for each one.
    ! Comment lines (starting with '!') are optional.
    !
    ! Gas units (must be same for all profiles)
    ! 0 => ppmv over dry air
    ! 1 => kg/kg over moist air
    ! 2 => ppmv over moist air
    !
    2
    !
    """

    for jj in range(jjmax):
        for ii in range(iimax):
            print("Creating profile data for each grid point within the WRF domain ..", end="\r")
            profile_file = f"prof-{profileCount:06}.dat"

            with open(profile_file, "w") as f:

                f.write(header_text)

                f.write(f"! --- Profile {profileCount} ---\n")

                # ===== SORT LEVELS (TOP → SURFACE) =====
                p_col = p_half[:, jj, ii]
                sort_idx = np.argsort(p_col)  # ascending

                p_sorted = p_col[sort_idx]
                T_sorted = temperature[:, jj, ii][sort_idx[:-1]]
                q_sorted = qv[:, jj, ii][sort_idx[:-1]]

                # --- Pressure ---
                f.write("!\n! Pressure levels (hPa)\n!\n")
                for val in p_sorted:
                    f.write(f"{val/100.0:.6f}\n")

                # --- Temperature ---
                f.write("!\n! Temperature profile (K)\n!\n")
                for val in T_sorted:
                    f.write(f"{val:.6f}\n")

                # --- Humidity ---
                f.write("!\n! Water vapour profile (kg/kg)\n!\n")
                for val in q_sorted:
                    f.write(f"{max(val,1e-5):.6f}\n")

                # --- Near surface ---
                f.write("!\n! Near-surface variables:\n")
                f.write("!  2m T (K)    2m q (ppmv)  10m wind u (m/s)  10m wind v (m/s)  wind fetch (m)\n!\n")
                f.write(f"   {t2m[jj,ii]:.4f}    {q2m[jj,ii]:.4f}     {u10[jj,ii]:.3f}             {v10[jj,ii]:.4f}            100000.\n!\n")

                # --- Skin ---
                f.write("! Skin variables:\n")
                f.write("!  Skin T (K)  Salinity   FASTEM parameters for land surfaces\n!\n")
                f.write(f"   {skinT[jj,ii]:.4f}    35.0       3.0 5.0 15.0 0.1 0.3\n!\n")

                # --- Surface ---
                f.write("! Surface type (0=land, 1=sea, 2=sea-ice) and water type (0=fresh, 1=ocean)\n!\n")
                f.write(f"   {int(landSeaMask[jj,ii])}         1\n!\n")

                # --- Elevation ---
                f.write("! Elevation (km), latitude and longitude (degrees)\n!\n")
                altitude = geopotential[jj, ii] / 9810.0
                f.write(f"   {altitude:.3f}    {lat[jj]:.3f}   {lon[ii]:.3f}\n!\n")

                # --- Angles ---
                f.write("! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n!\n")

                obs = get_observer_look(
                    satLon, satLat, satAltitude,
                    observationTime,
                    np.array([lon[ii]]),
                    np.array([lat[jj]]),
                    np.array([altitude])
                )

                satZen = 90 - obs[1][0]
                satAz = obs[0][0]

                sun = get_alt_az(observationTime, lon[ii], lat[jj])
                sunZen = sun[0] * 180 / pi
                sunAz = sun[1] * 180 / pi

                f.write(f"   {satZen:.3f}     {satAz:.3f}     {sunZen:.3f}     {sunAz:.3f}\n!\n")

                # --- Cloud ---
                f.write("! Cloud top pressure (hPa) and cloud fraction for simple cloud scheme\n!\n")
                f.write(f"   500.00    {cloudFraction[jj,ii]:.3f}\n")

                f.write("!\n! --- End of profile ---\n")

            os.rename(profile_file, os.path.join(dirName, profile_file))
            profileCount += 1
    print("Profile files have been created for RTTOV. Total profiles:", profileCount-1)
    application_shell.make_final_application_shell(rttovCoef, str(nlev+1), str(do_solar), satChannels, str(len(satChannels000)), rttov_install_path)

def make_netcdf():
    try:
        importlib.import_module("xarray")
    except:
        print("Warning: The Python module xarray is not installed.")
        print("Install it and run again.")
        print("Exiting ..")
        exit()
    import re
    import xarray as xr
    from glob import glob
    if os.path.exists(era5_surface_file000[0]) and os.path.exists(era5_level_file000[0]):
            era5_surface_file = era5_surface_file000[0]
            era5_level_file = era5_level_file000[0]
            print(era5_surface_file, "and", era5_level_file)
    else:
        print(f"Two ERA5 data files starting from era5data_surface_level and era5data_pressure_levels are missing")
        print("You can run again the application if you don't have the files")
        exit()
    # era5filexr = xr.open_dataset(era5_surface_file, engine='netcdf4', mode='r')
    era5filexr = xr.open_dataset(era5_surface_file, engine='cfgrib',
                                 backend_kwargs={"filter_by_keys": {"typeOfLevel": "surface"}, "indexpath": ""})
    allOutputs = glob(filePrefix+"_outputs/*")
    t2 = era5filexr.t2m
    # t2 = t2000.squeeze('valid_time')
    xlat  = era5filexr.latitude.to_numpy()
    xlong  = era5filexr.longitude.to_numpy()
    fillvalue = np.float32(-9999)
    jjmax = xlat.size
    iimax = xlong.size
    fillerVar = t2.to_numpy()
    fillerVar[:] = 0
    brightnessTemperature = xr.Dataset(
        coords={
            "latitude": (["latitude"], xlat),  # 2D latitude coordinates
            "longitude": (["longitude"], xlong),  # 2D longitude coordinates
        },
        attrs={
            "_FillValue": fillvalue,
            "title": "brightness_temperature",
            },
    )

    for myband in bandNames:
        brightnessTemperature[myband] = (["latitude", "longitude"], fillerVar.copy()) #copy is critical. Otherwise, all variables would be modified
        brightnessTemperature[myband].attrs["units"] = "K"
        brightnessTemperature[myband].attrs["long_name"] = "Calculated brightness temperature"
        brightnessTemperature[myband].attrs["_FillValue"] = fillvalue

    radiance = brightnessTemperature.copy(deep=True) #defining the variables coordinates, as well as filling them
    radiance.attrs["title"] = "radiance"
    overcastRadiance = brightnessTemperature.copy(deep=True)
    overcastRadiance.attrs["title"] = "overcast_radiances"
    transmittance = brightnessTemperature.copy(deep=True)
    transmittance.attrs["title"] = "transmittance"
    emissivity = brightnessTemperature.copy(deep=True)
    emissivity.attrs["title"] = "emissivities"
    for myband in bandNames: #modifying units and file titles
        radiance[myband].attrs["units"] = "mW/m2/sr/cm-1"
        radiance[myband].attrs["long_name"] = "Calculated radiance"
        overcastRadiance[myband].attrs["units"] = ""
        overcastRadiance[myband].attrs["long_name"] = "Calculated overcast radiances"
        transmittance[myband].attrs["units"] = ""
        transmittance[myband].attrs["long_name"] = "Calculated surface to space transmittance"
        emissivity[myband].attrs["units"] = ""
        emissivity[myband].attrs["long_name"] = "Calculated surface emissivities"
    print("Extracting the RTTOV outputs and storing them in arrays ..")
    for jj in range(jjmax): #latitudesTemperature profile (K)
        for ii in range(iimax): #longitude
            counter = jj*iimax + ii
            with open(allOutputs[counter], "r") as file:
                first_line = file.readline()
                if re.search("missing_value", first_line):
                    for band in bandNames:
                        brightnessTemperature[band][jj, ii] = fillvalue
                        radiance[band][jj, ii] = fillvalue
                        overcastRadiance[band][jj, ii] = fillvalue
                        transmittance[band][jj, ii] = fillvalue
                        emissivity[band][jj, ii] = fillvalue
                else:
                    found = False
                    fileExtract = []
                    lines = file.readlines()
                    for line in lines:
                        if "CHANNELS PROCESSED" in line:
                            found = True
                        if found:
                            fileExtract.append(line.strip())
                    for lineNumber, lineText in enumerate(fileExtract):
                        if "BRIGHTNESS TEMPERATURES" in lineText:
                            extractedLine = fileExtract[lineNumber + 1].strip()
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                brightnessTemperature[bandNames[chIndex]][jj, ii] = np.float32(simulationValue)
                            continue
                        if "CALCULATED RADIANCES" in lineText:
                            extractedLine = fileExtract[lineNumber + 1].strip()
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                radiance[bandNames[chIndex]][jj, ii] = np.float32(simulationValue)
                            continue
                        if "OVERCAST RADIANCES" in lineText:
                            extractedLine = fileExtract[lineNumber + 1].strip()
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                overcastRadiance[bandNames[chIndex]][jj, ii] = np.float32(simulationValue)
                            continue
                        if "SURFACE TO SPACE TRANSMITTANCE" in lineText:
                            extractedLine = fileExtract[lineNumber + 1].strip()
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                transmittance[bandNames[chIndex]][jj, ii] = np.float32(simulationValue)
                            continue
                        if "EMISSIVITIES" in lineText:
                            extractedLine = fileExtract[lineNumber + 1].strip()
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                emissivity[bandNames[chIndex]][jj, ii] = np.float32(simulationValue)

    print("Storing extracted values in NetCDF files ..")
    os.makedirs(postprocessingDir, exist_ok=True)
    print("Storing", brightnessTemperature.title, "into NetCDF")
    brightnessTemperature.to_netcdf(os.path.join(postprocessingDir, brightnessTemperature.title+".nc") )
    print("Storing", radiance.title, "into NetCDF")
    radiance.to_netcdf(os.path.join(postprocessingDir, radiance.title+".nc") )
    print("Storing", overcastRadiance.title, "into NetCDF")
    overcastRadiance.to_netcdf(os.path.join(postprocessingDir, overcastRadiance.title+".nc") )
    print("Storing", transmittance.title, "into NetCDF")
    transmittance.to_netcdf(os.path.join(postprocessingDir, transmittance.title+".nc") )
    print("Storing", emissivity.title, "into NetCDF")
    emissivity.to_netcdf(os.path.join(postprocessingDir, emissivity.title+".nc") )
    print("Brightness temperature, Radiance, Overcast radiance, Surface to space transmittance," + \
            "and emissivities have been stored in NetCDF files")

def plot_png():
    required_modules = ["matplotlib", "cartopy"]
    for module in required_modules:
        try:
            importlib.import_module(module)
        except:
            print("Warning: The Python module", module, "is not installed.")
            print("Install it and run again.")
            print("Exiting ..")
            exit()
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    for ncFile in glob(postprocessingDir+"/*.nc"):
        ds = xr.open_dataset(ncFile)
        lat = ds["latitude"]
        lon = ds["longitude"]
        subplot_parameters = {"projection": ccrs.PlateCarree()}
        for band in ds.data_vars.keys():
            print("Plotting", band, "of", ncFile)
            myvar = ds[band].copy(deep=True)
            fig, ax = plt.subplots( subplot_kw=subplot_parameters, figsize=(8, 6) )
            temp_plot = ax.pcolormesh(lon, lat, myvar, transform=ccrs.PlateCarree(), cmap="coolwarm")
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linestyle=":")
            gl = ax.gridlines(draw_labels=True, linestyle="-", alpha=1, color="black", linewidth=0.3)
            gl.xlocator = plt.MaxNLocator(8)  # Customize the number of longitude lines
            gl.ylocator = plt.MaxNLocator(6)  # Customize the number of latitude lines
            gl.top_labels = True  # Show tick marks on the top
            gl.right_labels = True  # Show tick marks on the right
            gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize x-axis tick labels
            gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize y-axis tick labels
            # Optional: Customize graticule labels format (degrees and direction)
            gl.xformatter = plt.FuncFormatter(lambda x, _: f"{x}°")
            gl.yformatter = plt.FuncFormatter(lambda y, _: f"{y}°")
            plt.colorbar(temp_plot, ax=ax, orientation="vertical", label=myvar.attrs["units"])
            ax.set_title(ds.attrs["title"] + " - " + band)
            plt.savefig(postprocessingDir + "/" + ds.attrs["title"]+ "_" + band + ".png")

def plot_rgb():
    required_modules = ["matplotlib", "cartopy"]
    for module in required_modules:
        try:
            importlib.import_module(module)
        except:
            print("Warning: The Python module", module, "is not installed.")
            print("Install it and run again.")
            print("Exiting ..")
            exit()
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import re
    redPol = namelist["postprocessing"]['RGB_plot_brightness_temperature']['Red']
    greenPol = namelist["postprocessing"]['RGB_plot_brightness_temperature']['Green']
    bluePol = namelist["postprocessing"]['RGB_plot_brightness_temperature']['Blue']
    redBands = set(re.findall(r'\b[a-zA-Z_]\w*\b', redPol))
    greenBands = set(re.findall(r'\b[a-zA-Z_]\w*\b', greenPol))
    blueBands = set(re.findall(r'\b[a-zA-Z_]\w*\b', bluePol))
    ncFile = postprocessingDir+"/brightness_temperature.nc"
    ds = xr.open_dataset(ncFile)
    for band in redBands:
        globals()[band] = ds[band].copy(deep=True)
    redPolXarray = eval(redPol)
    for band in greenBands:
        globals()[band] = ds[band].copy(deep=True)
    greenPolXarray = eval(greenPol)
    for band in blueBands:
        globals()[band] = ds[band].copy(deep=True)
    bluePolXarray = eval(bluePol)
    normalCoef = max([np.max(redPolXarray), np.max(greenPolXarray), np.max(bluePolXarray)])
    rgb_image = np.stack([redPolXarray/normalCoef, greenPolXarray/normalCoef, bluePolXarray/normalCoef], axis=-1)
    lat = ds["latitude"]
    lon = ds["longitude"]
    subplot_parameters = {"projection": ccrs.PlateCarree()}
    print("Plotting RGB image of brightness_temperature ..")
    fig, ax = plt.subplots( subplot_kw=subplot_parameters, figsize=(8, 6) )
    temp_plot = ax.pcolormesh(lon, lat, rgb_image, transform=ccrs.PlateCarree(), cmap="coolwarm")
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    gl = ax.gridlines(draw_labels=True, linestyle="-", alpha=1, color="black", linewidth=0.3)
    gl.xlocator = plt.MaxNLocator(8)  # Customize the number of longitude lines
    gl.ylocator = plt.MaxNLocator(6)  # Customize the number of latitude lines
    gl.top_labels = True  # Show tick marks on the top
    gl.right_labels = True  # Show tick marks on the right
    gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize x-axis tick labels
    gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize y-axis tick labels
    # Optional: Customize graticule labels format (degrees and direction)
    gl.xformatter = plt.FuncFormatter(lambda x, _: f"{x}°")
    gl.yformatter = plt.FuncFormatter(lambda y, _: f"{y}°")
    plt.colorbar(temp_plot, ax=ax, orientation="vertical", label=globals()[band].attrs["units"])
    plt.title("RGB Image from Brightness Temperatures")
    plt.axis("off")
    plt.savefig(postprocessingDir + "/" + "brightness_temperature_rgb.png")

# if __name__ == "__main__":
if not postprocessingEnabled:
    print("Postprocessing is disabled. Running to make the RTTOV profile files ..")
    make_inputdata()
else:
    postprocessing_directory_suffix = namelist["postprocessing"]['postprocessing_directory_suffix']
    postprocessingDir = filePrefix+"_"+postprocessing_directory_suffix
    rttov_output_dir = filePrefix+"_outputs"
    image_plot_enabled = namelist["postprocessing"]['image_plot_all_bands']
    rgb_plot_enabled = namelist["postprocessing"]['RGB_plot_brightness_temperature']['enabled']
    print("Postprocessing ...")
    if os.path.isdir(postprocessingDir): # NetCDF directory exists
        print(f"The postprocessing directory ({postprocessingDir}) already exists and may contain the NetCDF files.")
        print("Skipping making NetCDF files from RTTOV outputs.")
    else: # Making NetCDF files
        if not os.path.isdir(rttov_output_dir):
            print(f"Warning: The RTTOV outputs directory ({rttov_output_dir}) is not availabe.")
            print("Disable postprocessing and run to make the profile files.")
            print("Exiting ..")
            exit()
        print(f"Converting the RTTOV output within the ({rttov_output_dir}) directory to NetCDF files ..")
        make_netcdf()
    if image_plot_enabled:
        if os.path.isdir(postprocessingDir):
            print(f"Looking for the RTTOV NetCDF outputs in {postprocessingDir} ..")
            plot_png()
    if rgb_plot_enabled:
        if os.path.isdir(postprocessingDir):
            print(f"Looking for the RTTOV NetCDF outputs for brightness temperature in {postprocessingDir} ..")
            plot_rgb()