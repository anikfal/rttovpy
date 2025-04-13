## Code for downloading ERA5 data,and making RTTOV profile data based on the downloaded data,
## As well as generating the shellscript for running the RTTOV model using the profile data
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################
import importlib
required_modules = ["numpy", "yaml", "netCDF4", "pyorbital"]
for module in required_modules:
    try:
        importlib.import_module(module)
    except:
        print("Warning: The Python module", module, "is not installed.")
        print("Install it and run again.")
        print("Exiting ..")
        exit()

from modules import era5_download_manager, count_lines, application_shell
from modules.conversions import surface_humidity
import yaml, os
import netCDF4 as nc
from glob import glob
import numpy as np
from pyorbital.orbital import get_observer_look
from pyorbital.orbital import Orbital
from pyorbital.astronomy import get_alt_az
from datetime import datetime
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
angleEnable = namelist["satellite_information"]["user_defined_position"]['enabled']

with open('satellite_names.yaml', 'r') as yaml_file:
    satNameFile = yaml.safe_load(yaml_file)

filePrefix = namelist["area_of_simulation"]["domain_name"]
dirName = filePrefix+"_profiles/"
era5_surface_file000 = glob("era5data_surface_level_*" + filePrefix + ".nc")#[0]
era5_level_file000 = glob("era5data_pressure_levels_*" + filePrefix + ".nc")#[0]
rttovCoef = namelist["rttov_coefficient_file_path"]
postprocessingEnabled = namelist["postprocessing"]['enabled']

def make_inputdata():
    if os.path.exists(dirName) and os.path.isdir(dirName):
        print("Using the previously downloaded data in", dirName, "to make the final shellscript application (run_era5_example_fwd.sh)")
        profilesList000 = os.listdir(dirName)
        profilesList = [var for var in profilesList000 if (var.startswith("prof-") and var.endswith(".dat"))]
        pressureLevelsSize = count_lines.count_lines_between(dirName+profilesList[0], "! Pressure levels (hPa)", "! Temperature profile (K)") - 2
        application_shell.make_final_application_shell(rttovCoef, str(pressureLevelsSize), satChannels, rttov_install_path)
        exit()

    try:
        if os.path.exists(era5_surface_file000[0]) and os.path.exists(era5_level_file000[0]):
            era5_surface_file = era5_surface_file000[0]
            era5_level_file = era5_level_file000[0]
            print(f"  The files {era5_surface_file} and {era5_level_file} already exist and will be used to make profile data.")
            print(f"  So no new ERA5 data will be downloaded.")
    except:
            era5_download_manager.main_dm()
            era5_surface_file = glob("era5data_surface_level_*" + filePrefix + ".nc")[0]
            era5_level_file = glob("era5data_pressure_levels_*" + filePrefix + ".nc")[0]

    try:
        os.makedirs(dirName)
        print(f"Directory {dirName} has been created to store profile datafiles.")
    except Exception as error:
        print(f"An error occurred while creating {dirName}: {error}")

    year = namelist["time_of_simulation"]["year"]
    month = namelist["time_of_simulation"]["month"]
    day = namelist["time_of_simulation"]["day"]
    hour = namelist["time_of_simulation"]["hour"]

    orb = Orbital(satNameFile[satIndex])
    observationTime = datetime(year, month, day, hour)
    # satPositions= orb.get_lonlatalt(observationTime) #Get longitude, latitude and altitude of the satellite
    satPositions = (136.85902196460546, -53.70781534686423, 715.6113205704698)
    satAltitude = satPositions[2]
    if(angleEnable):
        print("Simulation with user defined satellite position:")
        satLat = namelist["satellite_information"]["user_defined_position"]['sat_latitude']
        satLon = namelist["satellite_information"]["user_defined_position"]['sat_longitude']
    else:
        satLat = satPositions[1]
        satLon = satPositions[0]

    ncFile_surface = nc.Dataset(era5_surface_file)
    t2m = np.round(ncFile_surface["t2m"], 4)
    d2m = np.round(ncFile_surface["d2m"], 4)
    p2m = np.round(ncFile_surface["sp"], 4)
    u10 = np.round(ncFile_surface["u10"], 4)
    v10 = np.round(ncFile_surface["v10"], 4)
    skinT = np.round(ncFile_surface["skt"], 4)
    landSeaMask = np.round(ncFile_surface["lsm"], 0)
    lakeCover = ncFile_surface["cl"]
    soilType = ncFile_surface["slt"]
    geopotential = np.round(ncFile_surface["z"], 4)
    # cloudBase = ncFile_surface["cbh"]
    cloudFraction = ncFile_surface["tcc"]
    q2m = np.round(surface_humidity(t2m[:], p2m[:]), 4)

    lat = ncFile_surface['latitude']
    lon = ncFile_surface['longitude']

    ncFile_level = nc.Dataset(era5_level_file)
    temperature000 = np.round(ncFile_level["t"], 4)
    qv = np.round(ncFile_level["q"], 5)
    varShape = qv.shape
    tempLevel = np.round(ncFile_level["t"], 4)
    levelVar = list(ncFile_level.dimensions)[1]
    pressureLevels = ncFile_level[levelVar]

    tt = 0
    profileCount = 1
    jjmax = varShape[2]
    iimax = varShape[3]
    for jj in range(jjmax): #latitudesTemperature profile (K)
        for ii in range(iimax): #longitude
    # for jj in range(5, 10): #latitudesTemperature profile (K)
    #     for ii in range(iimax): #longitude
            jjcount = jj+1
            iicount = ii+1
            print("Creating profile data for the grid point jj:", jjcount, "ii:", iicount)
            profileIndexNaming = f"j and i = {jjcount}/{jjmax} and {iicount}/{iimax}"
            profile_file = f"prof-{profileCount:06}.dat"
            profileCount = profileCount + 1
            with open(profile_file, "w") as file_writer:
                header_lines = [
                    f"! jj, ii, {jjcount}, {iicount}\n",
                    f"! jjSize, iiSize, {jjmax}, {iimax}\n",
                    "! Specify input profiles for example_fwd.F90.\n",
                    "! Multiple profiles may be described: follow the same format for each one.\n",
                    "! Comment lines (starting with '!') are optional.\n",
                    "!\n",
                    "! Gas units (must be same for all profiles)\n",
                    "! 0 => ppmv over dry air\n",
                    "! 1 => kg/kg over moist air\n",
                    "! 2 => ppmv over moist air\n",
                    "!\n",
                    " 1\n",
                    "!\n",
                    ]
                file_writer.writelines(header_lines)
            file_append = open(profile_file, "a")
            file_append.write("! --- Start of profile ---" + "\n")
            file_append.write("! --- Grid point with " + profileIndexNaming + " ---" + "\n")
            subHead = [
                    "!\n",
                    "! Pressure levels (hPa)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            levelRange = list(range(varShape[1]))[::-1]
            for level in levelRange:
                file_append.write(str(pressureLevels[level])+'\n')
            subHead = [
                    "!\n",
                    "! Temperature profile (K)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            for level in levelRange:
                pointValue = tempLevel[tt,level,jj,ii]
                file_append.write(str(pointValue)+'\n')
            subHead = [
                    "!\n",
                    "! Water vapour profile (kg/kg)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            for level in levelRange:
                pointValue = qv[0,level,jj,ii]
                if pointValue < 0.00001:
                    pointValue = np.float64(0.00001)
                file_append.write(np.array2string(pointValue, formatter={'float_kind':lambda x: f"{x:.{5}f}"})+'\n')
            subHead = [
                    "!\n",
                    "! Near-surface variables:\n"
                    "!  2m T (K)    2m q (kg/kg) 2m p (hPa) 10m wind u (m/s)  10m wind v (m/s)  wind fetch (m)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            near_surface_vars = [t2m[tt, jj, ii], q2m[tt, jj, ii], p2m[tt, jj, ii]/100,\
                                u10[tt, jj, ii], v10[tt, jj, ii], 100000]
            near_surface_vars_2line = ' '.join(map(str, near_surface_vars))
            file_append.write(near_surface_vars_2line)
            subHead = [
                    "\n!\n",
                    "! Skin variables:\n"
                    "!  Skin T (K)  Salinity   FASTEM parameters for land surfaces\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            skinVars = [skinT[tt, jj, ii], 35.0, 3.0, 5.0, 15.0, 0.1, 0.3]
            skinVars_2line = ' '.join(map(str, skinVars))
            file_append.write(skinVars_2line)
            subHead = [
                    "\n!\n",
                    "! Surface type (0=land, 1=sea, 2=sea-ice) and water type (0=fresh, 1=ocean)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            surfaceType = [int(landSeaMask[tt, jj, ii]), 1]
            surfaceType_2line = ' '.join(map(str, surfaceType))
            file_append.write(surfaceType_2line)
            subHead = [
                    "\n!\n",
                    "! Elevation (km), latitude and longitude (degrees)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            observerAltitude = np.round(geopotential[tt, jj, ii]/9810, 4)
            elevation = [observerAltitude, lat[jj], lon[ii]]
            elevation_2line = ' '.join(map(str, elevation))
            file_append.write(elevation_2line)
            subHead = [
                    "\n!\n",
                    "! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            userdefSatPos = get_observer_look(satLon, satLat, satAltitude, observationTime, np.array([lon[ii]]), np.array([lat[jj]]), np.array([observerAltitude]))
            satAzimuth = userdefSatPos[0]
            satZenith = 90 - userdefSatPos[1]
            sunPositions = get_alt_az(observationTime, lon[ii], lat[jj])
            sunZenith = sunPositions[0] * 180/pi
            sunAzimuth = sunPositions[1] * 180/pi
            satsunAngles = np.round([satZenith[0], satAzimuth[0], sunZenith, sunAzimuth], 2)
            satsunAngles_2line = ' '.join(map(str, satsunAngles))
            file_append.write(satsunAngles_2line)
            subHead = [
                    "\n!\n",
                    "! Cloud top pressure (hPa) and cloud fraction for simple cloud scheme\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            cloudinfo = [500, cloudFraction[tt, jj, ii]]
            cloudinfo_2line = ' '.join(map(str, cloudinfo))
            file_append.write(cloudinfo_2line)
            subHead = [
                    "\n!\n",
                    "! --- End of profile " + profileIndexNaming + " ---" + "\n"
                    "!\n",
                ]
            for line in subHead:
                file_append.write(line)
            destination_path = os.path.join(dirName, os.path.basename(profile_file))
            try:
                os.rename(profile_file, destination_path)
            except Exception as error:
                print(f"An error occurred: {error}")
                print("Exiting ..")
                exit()
            # exit()
            file_append.close()

    application_shell.make_final_application_shell(rttovCoef, str(varShape[1]), satChannels, rttov_install_path)

import xarray as xr
from glob import glob
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
    era5filexr = xr.open_dataset(era5_surface_file, engine='netcdf4', mode='r')
    allOutputs = glob(filePrefix+"_outputs/*")
    t2000 = era5filexr.t2m
    t2 = t2000.squeeze('valid_time')
    xlat  = era5filexr.latitude.to_numpy()
    xlong  = era5filexr.longitude.to_numpy()
    fillvalue = np.float32(-9999)
    jjmax = xlat.size
    iimax = xlong.size
    fillerVar = t2000.to_numpy()[0,:,:]
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
            counter = jj*jjmax + ii
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

if __name__ == "__main__":
    if not postprocessingEnabled:
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