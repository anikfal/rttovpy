## Code for downloading ERA5 data,and making RTTOV profile data based on the downloaded data,
## As well as generating the shellscript for running the RTTOV model using the profile data
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################

from modules import count_lines, application_shell
import yaml, os
import netCDF4 as nc
import numpy as np
from pyorbital.orbital import get_observer_look
from pyorbital.orbital import Orbital
from pyorbital.astronomy import get_alt_az
from datetime import datetime, timedelta
from math import pi
from wrf import getvar, disable_xarray
disable_xarray()


with open('namelist_wrf.yaml', 'r') as yaml_file:
    namelist = yaml.safe_load(yaml_file)
satIndex = namelist["satellite_information"]["sat_name_index"]
satChannels000 = namelist["satellite_information"]["sat_channel_list"]
satChannels = " ".join(str(item_value) for item_value in satChannels000)
rttov_install_path = namelist["rttov_installation_path"]
if rttov_install_path.endswith('/'):
    rttov_install_path = rttov_install_path[:-1]
if satIndex > 76:
    print("Warning: Satellite name index must be between 1 to 76. Please look inside <satellite_names.yaml>.")
    print("exiting ..")
    exit()
angleEnable = namelist["satellite_information"]["user_defined_position"]['enable']

with open('satellite_names.yaml', 'r') as yaml_file:
    satNameFile = yaml.safe_load(yaml_file)

wrfFilePath = namelist["wrf_file_path"]
wrfFileName = os.path.basename(wrfFilePath)

if not os.path.exists(wrfFilePath):
    print("Warning:", wrfFilePath, "is not a valid file path.")
    print("Exiting ..")

dirnameSuffix = namelist["profiles_direname_suffix"]
dirName = wrfFileName+dirnameSuffix+"/"
rttovCoef = namelist["rttov_coefficient_file_path"]

if os.path.exists(dirName) and os.path.isdir(dirName):
    print("- Using the previously downloaded data in", dirName, "to make the final RTTOV shell application")
    profilesList000 = os.listdir(dirName)
    profilesList = [var for var in profilesList000 if (var.startswith("prof-") and var.endswith(".dat"))]
    pressureLevelsSize = count_lines.count_lines_between(dirName+profilesList[0], "! Pressure levels (hPa)", "! Temperature profile (K)") - 2
    application_shell.make_final_application_shell(rttovCoef, str(pressureLevelsSize), satChannels, rttov_install_path)
    print("- The file 'run_era5_example_fwd.sh' has been made")
    exit()

try:
    pass
    os.makedirs(dirName)
    print(f"Directory {dirName} has been created to store profile datafiles.")
except Exception as error:
    print(f"An error occurred while creating {dirName}: {error}")

wrffile = nc.Dataset(wrfFilePath)
minuteArr = wrffile.variables["XTIME"]
startDate000 = wrffile.START_DATE
startYear = int(startDate000.split("-")[0])
startMonth = int(startDate000.split("-")[1])
startDay = int(startDate000.split("-")[2][:2])
startHour = int(startDate000.split("-")[2][3:5])
startDate = datetime(startYear, startMonth, startDay, startHour)
endDate = startDate + timedelta(minutes=int(minuteArr[-1]))

year = namelist["time_of_simulation"]["year"]
month = namelist["time_of_simulation"]["month"]
day = namelist["time_of_simulation"]["day"]
hour = namelist["time_of_simulation"]["hour"]
observationTime = datetime(year, month, day, hour)
if not startDate <= observationTime <= endDate:
    print("!!!!!!!!!!!!!!!!!!!")
    print(f"Warning: {observationTime} is not among the time-slots of {wrfFileName} (must be between {startDate} and {endDate}).")
    print("Exiting ..")
    exit()

observationIndex000 = (observationTime - startDate)/ (60*int(minuteArr[1]))
observationIndex = int(observationIndex000.total_seconds()) - 1
# orb = Orbital(satNameFile[satIndex])
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

t2m = getvar(wrffile, "T2", timeidx=observationIndex)
d2m = getvar(wrffile, "td2", timeidx=observationIndex)+273.15
qv = getvar(wrffile, "QVAPOR", timeidx=observationIndex)
p = getvar(wrffile, "p", timeidx=observationIndex)
p2m = p[0,:,:]
u10 = getvar(wrffile, "U10", timeidx=observationIndex)
v10 = getvar(wrffile, "V10", timeidx=observationIndex)
skinT = getvar(wrffile, "TSK", timeidx=observationIndex)
lakeCover = getvar(wrffile, "LAKEMASK", timeidx=observationIndex)
modelheight = getvar(wrffile, "HGT", timeidx=observationIndex) #altitude
cloudFraction000 = getvar(wrffile, "cloudfrac", timeidx=observationIndex)
cloudFraction = np.maximum.reduce([cloudFraction000[0,:,:], cloudFraction000[1,:,:], cloudFraction000[2,:,:]]) #maximum values
q2m = qv[0,:,:]
lat = getvar(wrffile, "lat", timeidx=observationIndex)
lon = getvar(wrffile, "lon", timeidx=observationIndex)
tempLevel = getvar(wrffile, "tk", timeidx=observationIndex)

varShape = qv.shape
pressureLevels = p/100 #hpa

profileCount = 1
jjmax = varShape[1]
iimax = varShape[2]
for jj in range(jjmax): #latitudesTemperature profile (K)
    for ii in range(iimax): #longitude
# for jj in range(3): #latitudesTemperature profile (K)
#     for ii in range(2): #longitude
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
        levelRange = list(range(varShape[0]))[::-1]
        for level in levelRange:
            file_append.write(str(pressureLevels[level,jj,ii])+'\n')

        subHead = [
                "!\n",
                "! Temperature profile (K)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        for level in levelRange:
            # pointValue = tempLevel[tt,level,jj,ii]
            pointValue = tempLevel[level,jj,ii]
            file_append.write(str(pointValue)+'\n')
        
        subHead = [
                "!\n",
                "! Water vapour profile (kg/kg)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        for level in levelRange:
            pointValue = qv[level,jj,ii]
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
        near_surface_vars = [t2m[jj, ii], q2m[jj, ii], p2m[jj, ii]/100,\
                              u10[jj, ii], v10[jj, ii], 100000]
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
        skinVars = [skinT[jj, ii], 35.0, 3.0, 5.0, 15.0, 0.1, 0.3]
        skinVars_2line = ' '.join(map(str, skinVars))
        file_append.write(skinVars_2line)

        subHead = [
                "\n!\n",
                "! Surface type (0=land, 1=sea, 2=sea-ice) and water type (0=fresh, 1=ocean)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        surfaceType = [int(lakeCover[jj, ii]), 1]
        surfaceType_2line = ' '.join(map(str, surfaceType))
        file_append.write(surfaceType_2line)

        subHead = [
                "\n!\n",
                "! Elevation (km), latitude and longitude (degrees)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        # observerAltitude = np.round(modelheight[jj, ii]/9810, 4)
        observerAltitude = modelheight[jj, ii]/1000
        elevation = [observerAltitude, lat[jj, ii], lon[jj, ii]]
        elevation_2line = ' '.join(map(str, elevation))
        file_append.write(elevation_2line)

        subHead = [
                "\n!\n",
                "! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        userdefSatPos = get_observer_look(satLon, satLat, satAltitude, observationTime, lon[jj, ii], lat[jj, ii], observerAltitude)
        satAzimuth = userdefSatPos[0]
        satZenith = 90 - userdefSatPos[1]
        sunPositions = get_alt_az(observationTime, lon[jj, ii], lat[jj, ii])
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
        cloudinfo = [500, cloudFraction[jj, ii]]
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

        file_append.close()

print("")
print("==================================================================")
print("Making the shellscript application for the RTTOV forward model ...")
application_shell.make_final_application_shell(rttovCoef, str(varShape[0]), satChannels, rttov_install_path)
print("The file run_era5_example_fwd.sh has been made successfully.")