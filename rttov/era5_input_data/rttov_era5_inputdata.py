from modules import era5_download_manager, count_lines, application_shell
from modules.conversions import surface_humidity
import yaml
import netCDF4 as nc
from glob import glob
import numpy as np
from pyorbital.orbital import get_observer_look
from pyorbital.orbital import Orbital
from pyorbital.astronomy import get_alt_az
from datetime import datetime
from math import pi
import os, sys

with open('input.yaml', 'r') as yaml_file:
    inputFile = yaml.safe_load(yaml_file)
satIndex = inputFile["satellite_information"]["sat_name_index"]
satChannels000 = inputFile["satellite_information"]["sat_channel_list"]
satChannels = " ".join(str(item_value) for item_value in satChannels000)
rttov_install_path = inputFile["rttov_installation_path"]
if rttov_install_path.endswith('/'):
    rttov_install_path = rttov_install_path[:-1]
if satIndex > 76:
    print("Warning: Satellite name index must be between 1 to 76. Please look inside <satellite_names.yaml>.")
    print("exiting ..")
    exit()
angleEnable = inputFile["satellite_information"]["user_defined_position"]['enable']

with open('satellite_names.yaml', 'r') as yaml_file:
    satNameFile = yaml.safe_load(yaml_file)
# print("sat index:", satIndex, satNameFile[satIndex] )

if(angleEnable):
    print("Simulation with user defined satellite position:")
    satLat = inputFile["satellite_information"]["user_defined_position"]['sat_latitude']
    satLon = inputFile["satellite_information"]["user_defined_position"]['sat_longitude']

filePrefix = inputFile["area_of_simulation"]["domain_name"]
dirName = filePrefix+"_profiles/"
era5_surface_file = glob("era5data_surface_level_*" + filePrefix + ".nc")[0]
era5_level_file = glob("era5data_pressure_levels_*" + filePrefix + ".nc")[0]
rttovCoef = inputFile["rttov_coefficient_file_path"]





if os.path.exists(dirName) and os.path.isdir(dirName):
    print("Using the previously downloaded data in", dirName, "to make the final shellscript application (run_era5_example_fwd.sh)")
    varShape = count_lines.count_lines_between(dirName+"/prof-000001.dat", "! Pressure levels (hPa)", "! Temperature profile (K)") - 2
    application_shell.make_final_application_shell(rttovCoef, str(varShape), satChannels, rttov_install_path)
    exit()
if os.path.exists(era5_surface_file) and os.path.exists(era5_level_file):
    print(f"  The files {era5_surface_file} and {era5_level_file} already exist and will be used to make profile data.")
    print(f"  So no new ERA5 data will be downloaded.")
else:
    era5_download_manager.main_dm()




# with open('input.yaml', 'r') as yaml_file:
#     input_data = yaml.safe_load(yaml_file)
year = inputFile["time_of_simulation"]["year"]
month = inputFile["time_of_simulation"]["month"]
day = inputFile["time_of_simulation"]["day"]
hour = inputFile["time_of_simulation"]["hour"]

orb = Orbital(satNameFile[satIndex])
observationTime = datetime(year, month, day, hour)
satPositions= orb.get_lonlatalt(observationTime) #Get longitude, latitude and altitude of the satellite
# satPositions = (136.85902196460546, -53.70781534686423, 715.6113205704698)
satAlt = satPositions[2]

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


application_shell.make_final_application_shell()


# try:
#     os.makedirs(dirName)
#     print(f"Directory {dirName} has been created to store profile datafiles.")
# except FileExistsError:
#     print(f"Warning: directory {dirName} to store profile datafiles, already exist.")
#     print("Exiting ..")
#     exit()
# except Exception as error:
#     print(f"An error occurred: {error}")


tt = 0
profileCount = 1
jjmax = varShape[2]
iimax = varShape[3]
for jj in range(jjmax): #latitudesTemperature profile (K)
    for ii in range(iimax): #longitude
# for jj in range(5): #latitudesTemperature profile (K)
#     for ii in range(5): #longitude
        print("  jj:", jj, "ii:", ii)
        profileIndexNaming = f"j and i = {jj}/{jjmax} and {ii}/{iimax}"
        profile_file = f"prof-{profileCount:06}.dat"
        profileCount = profileCount + 1
        with open(profile_file, "w") as file_writer:
            header_lines = [
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
        userdefSatPos = get_observer_look(satLon, satLat, satAlt, observationTime, lon[ii], lat[jj], observerAltitude)
        satAzimuth = userdefSatPos[0]
        satZenith = 90 - userdefSatPos[1]
        sunPositions = get_alt_az(observationTime, lon[ii], lat[jj])
        sunZenith = sunPositions[0] * 180/pi
        sunAzimuth = sunPositions[1] * 180/pi
        satsunAngles = np.round([satZenith, satAzimuth, sunZenith, sunAzimuth], 2)
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