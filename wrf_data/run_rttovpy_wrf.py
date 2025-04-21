## Module for extracting WRF model data,and making RTTOV profile data based on the extracted data,
## As well as generating the shellscript for running the RTTOV model using the profile data
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################

import importlib
required_modules = ["numpy", "yaml", "netCDF4", "pyorbital", "wrf"]
for module in required_modules:
    try:
        importlib.import_module(module)
    except:
        print("Warning: The Python module", module, "is not installed.")
        print("Install it and run again.")
        print("Exiting ..")
        exit()

import yaml
import os
import netCDF4 as nc
from datetime import datetime, timedelta
import numpy as np
with open('namelist_wrf.yaml', 'r') as yaml_file:
    namelist = yaml.safe_load(yaml_file)
postprocessingEnabled = namelist["postprocessing"]['enabled']
dust = namelist["wrfchem_dust_profiles"]['enabled']
wrfFilePath = namelist["wrf_file_path"]
wrfFileName = os.path.basename(wrfFilePath)
wrffile = nc.Dataset(wrfFilePath)
minuteArr = wrffile.variables["XTIME"]
startDate000 = wrffile.START_DATE
startYear = int(startDate000.split("-")[0])
startMonth = int(startDate000.split("-")[1])
startDay = int(startDate000.split("-")[2][:2])
startHour = int(startDate000.split("-")[2][3:5])
startDate = datetime(startYear, startMonth, startDay, startHour)
endDate = startDate + timedelta(minutes=int(minuteArr[-1]))
dirnameSuffix = namelist["rttov_inputdata_directory_suffix"]
dirName = wrfFileName+"_"+dirnameSuffix+"/"

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
satChannels000 = namelist["satellite_information"]["sat_channel_list"]
satChannels = " ".join(str(item_value) for item_value in satChannels000)
bandNames = namelist["satellite_information"]["sat_channel_names"]
if len(satChannels000) != len(bandNames):
    print("Warning: lengths of sat_channel_list and sat_channel_names are not equal.")
    print("Exiting ..")
    exit()

def make_inputdata():
    from modules import count_lines, application_shell
    from pyorbital.orbital import get_observer_look
    from pyorbital.orbital import Orbital
    from pyorbital.astronomy import get_alt_az
    from math import pi
    from wrf import getvar, disable_xarray
    disable_xarray()

    satIndex = namelist["satellite_information"]["sat_name_index"]
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

    if not os.path.exists(wrfFilePath):
        print("Warning:", wrfFilePath, "is not a valid file path.")
        print("Exiting ..")
        exit()

    rttovCoef = namelist["rttov_coefficient_file_path"]

    # if os.path.exists(dirName) and os.path.isdir(dirName):
    if os.path.isdir(dirName) and bool(os.listdir(dirName)):
        print("- Using the previously downloaded data in", dirName, "to make the final RTTOV shell application")
        profilesList000 = os.listdir(dirName)
        profilesList = [var for var in profilesList000 if (var.startswith("prof-") and var.endswith(".dat"))]
        pressureLevelsSize = count_lines.count_lines_between(dirName+profilesList[0], "! Pressure levels (hPa)", "! Temperature profile (K)") - 2
        application_shell.make_final_application_shell(rttovCoef, str(pressureLevelsSize), satChannels, rttov_install_path)
        print("- The file 'run_wrf_example_fwd.sh' has been made")
        exit()

    try:
        pass
        os.makedirs(dirName)
        print(f"Directory {dirName} has been created to store profile datafiles.")
    except Exception as error:
        print(f"An error occurred while creating {dirName}: {error}")

    orb = Orbital(satNameFile[satIndex])
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
            userdefSatPos = get_observer_look(satLon, satLat, satAltitude, observationTime, np.array([lon[jj, ii]], dtype=np.float32), np.array([lat[jj, ii]], dtype=np.float32), np.array([observerAltitude], dtype=np.float32))
            satAzimuth = userdefSatPos[0]
            satZenith = 90 - userdefSatPos[1]
            sunPositions = get_alt_az(observationTime, lon[jj, ii], lat[jj, ii])
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

    if not dust:
        print("")
        print("==================================================================")
        print("Making the shellscript application for the RTTOV forward model ...")
        application_shell.make_final_application_shell(rttovCoef, str(varShape[0]), satChannels, rttov_install_path)
        print("The file run_wrf_example_fwd.sh has been made successfully.")
    else:
        aerosol_coefficient_file_path = namelist["wrfchem_dust_profiles"]['aerosol_coefficient_file_path']
        if not os.path.exists(aerosol_coefficient_file_path):
            print("Warning:", aerosol_coefficient_file_path, "is not a valid file path.")
            print("Exiting ..")
            exit()
        make_dust_profile()
        print("")
        print("==================================================================")
        print("Making the shellscript application for the RTTOV forward model ...")
        application_shell.make_final_dust_application_shell(rttovCoef, str(varShape[0]), satChannels, rttov_install_path, aerosol_coefficient_file_path)
        print("The file run_wrf_example_fwd.sh has been made successfully.")

def make_dust_profile():
    from modules import count_lines, application_shell
    from math import pi
    from wrf import getvar, disable_xarray
    disable_xarray()

    try:
        mineral_nuc = getvar(wrffile, "DUST_1", timeidx=observationIndex) * 1e-9
        mineral_acc = (getvar(wrffile, "DUST_2", timeidx=observationIndex) + getvar(wrffile, "DUST_3", timeidx=observationIndex) + getvar(wrffile, "DUST_4", timeidx=observationIndex)) * 1e-9
        mineral_coa = getvar(wrffile, "DUST_5", timeidx=observationIndex) * 1e-9
        p = getvar(wrffile, "p", timeidx=observationIndex)

    except:
        print(f"Warning: There is no dust variabls (DUST_1, DUST_2, ..., DUST_5) within {wrfFilePath}")
        print("This file is probably not a WRF/Chem output file.")
        print("Exiting ..")
        exit()

    varShape = mineral_nuc.shape

    profileCount = 1
    jjmax = varShape[1]
    iimax = varShape[2]
    fillerNull = "0.000"
    for jj in range(jjmax): #latitudesTemperature profile (K)
        for ii in range(iimax): #longitude
    # for jj in range(3): #latitudesTemperature profile (K)
    #     for ii in range(2): #longitude
            jjcount = jj+1
            iicount = ii+1
            print("Creating profile data for the grid point jj:", jjcount, "ii:", iicount)
            profileIndexNaming = f"j and i = {jjcount}/{jjmax} and {iicount}/{iimax}"
            profile_file = f"aer_prof-{profileCount:06}.dat"
            profileCount = profileCount + 1

            with open(profile_file, "w") as file_writer:
                header_lines = [
                    f"! jj, ii, {jjcount}, {iicount}\n",
                    f"! jjSize, iiSize, {jjmax}, {iimax}\n",
                    "! Specify input profiles for example_fwd.F90.\n",
                    "!\n", 
                    "! NB This file must contain data for the same number of profiles as prof.dat\n"
                    "!\n"
                    "! For each profile, the first value is an integer. If the value is 1-10\n"
                    "! this indicates that a climatological aerosol profile generated by\n"
                    "! rttov_aer_clim_prof should be used:\n"
                    "!\n"
                    "!         1  -->Continental clean\n"
                    "!         2  -->Continental average\n"
                    "!         3  -->Continental polluted\n"
                    "!         4  -->Urban\n"
                    "!         5  -->Desert\n"
                    "!         6  -->Maritime clean\n"
                    "!         7  -->Maritime polluted\n"
                    "!         8  -->Maritime tropical\n"
                    "!         9  -->Arctic\n"
                    "!         10 -->Antarctic\n"
                    "!\n"
                    "! Any other value implies that input profiles for all 13 aerosol types\n"
                    "! are provided here.\n"
                    "!\n"
                    "! Flag to indicate aerosol units (T => kg/kg; F => number density cm^-3):\n"
                    "!\n"
                    "  T\n"
                    "!\n"
                    ]
                file_writer.writelines(header_lines)
                file_append = open(profile_file, "a")
                file_append.write("! --- Start of profile ---" + "\n")
                file_append.write("! --- Grid point with " + profileIndexNaming + " ---" + "\n")
                file_append.write("! Supply the number density profiles here" + "\n")
                file_append.write("  0" + "\n")
                file_append.write("!" + "\n")
                file_append.write("! Dust concentration profiles (kg/kg) for each aerosol particle type (1-13) for each layer" + "\n")
                file_append.write("!" + "\n")
                file_append.write("!      INSO       WASO       SOOT       SSAM       SSCM       MINM       MIAM       MICM       MITR       SUSO       VOLA       VAPO       ASDU" + "\n")
                file_append.write("!" + "\n")

            levelRange = list(range(varShape[0]))[::-1]
            for level in levelRange:
                # row_value = ["      0.000", fillerNull, fillerNull, fillerNull, fillerNull, f"{mineral_nuc[level,jj,ii]:.3f}",
                            #  f"{mineral_nuc[level,jj,ii]:.3f}", f"{mineral_nuc[level,jj,ii]:.3f}", fillerNull, fillerNull, fillerNull, fillerNull, fillerNull]
                row_value = ["      0.000", fillerNull, fillerNull, fillerNull, fillerNull, str(mineral_nuc[level,jj,ii]),
                             str(mineral_acc[level,jj,ii]), str(mineral_coa[level,jj,ii]), fillerNull, fillerNull, fillerNull, fillerNull, fillerNull]
                row_value_str = "      ".join(row_value)
                file_append.write(row_value_str+'\n')
                destination_path = os.path.join(dirName, os.path.basename(profile_file))
            try:
                os.rename(profile_file, destination_path)
            except Exception as error:
                print(f"An error occurred: {error}")
                print("Exiting ..")
                exit()

            file_append.close()

import xarray as xr
from glob import glob
outputDirnameSuffix = namelist["rttov_outputdata_directory_suffix"]
basedir = os.path.basename(wrfFilePath)
outputDirPath = basedir+"_"+outputDirnameSuffix
def make_netcdf():
    import re
    wrffilexr = xr.open_dataset(wrfFilePath, engine='netcdf4', mode='r')
#    outputDirnameSuffix = namelist["rttov_outputdata_directory_suffix"]
#    basedir = os.path.basename(wrfFilePath)
#    outputDirPath = basedir+"_"+outputDirnameSuffix
    t2000 = wrffilexr.T2
    t2 = t2000.isel(Time=observationIndex).squeeze()
    xlat  = wrffilexr.XLAT.to_numpy()[0,:,:]
    xlong = wrffilexr.XLONG.to_numpy()[0,:,:]

    fillvalue = np.float32(-9999)
    varshape = t2.shape
    jjmax = varshape[0]
    iimax = varshape[1]
    allOutputs = glob(outputDirPath+"/output*")

    fillerVar = t2000.to_numpy()[0,:,:]
    fillerVar[:] = 0
    
    brightnessTemperature = xr.Dataset(
        coords={
            "lat": (["south_north", "west_east"], xlat),  # 2D latitude coordinates
            "lon": (["south_north", "west_east"], xlong),  # 2D longitude coordinates
        },
        attrs={
            "_FillValue": fillvalue,
            "title": "brightness_temperature",
            "projection": wrffilexr.attrs["MAP_PROJ_CHAR"],
            "central_latitude": wrffilexr.attrs["CEN_LAT"],
            "central_longitude": wrffilexr.attrs["CEN_LON"],
            "standard_parallels": (wrffilexr.attrs["TRUELAT1"], wrffilexr.attrs["TRUELAT2"]),
            },
    )
    for myband in bandNames:
        brightnessTemperature[myband] = (["south_north", "west_east"], fillerVar.copy()) #copy is critical. Otherwise, all variables would be modified
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
                            # print(extractedLine)
                            for chIndex, simulationValue in enumerate(extractedLine.split()):
                                # print("np.float32(simulationValue)", np.float32(simulationValue))
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
        lat = ds["lat"]
        lon = ds["lon"]
        projection = ds.attrs["projection"]
        if "lambert" in projection.lower():
            subplot_parameters = {"projection": ccrs.LambertConformal(central_longitude=ds.attrs["central_longitude"], 
                        central_latitude=ds.attrs["central_latitude"], standard_parallels=ds.attrs["standard_parallels"],)}
        elif "mercator" in projection.lower():
            subplot_parameters = {"projection": ccrs.Mercator( central_longitude=ds.attrs["central_longitude"], 
                        min_latitude=lat.min(), max_latitude=lat.max() )}
        elif "polar" in projection.lower():
            if lat.mean() > 0:
                subplot_parameters = {"projection": ccrs.NorthPolarStereo( central_longitude=ds.attrs["central_longitude"] )}
            else:
                subplot_parameters = {"projection": ccrs.SouthPolarStereo( central_longitude=ds.attrs["central_longitude"] )}
        else: #consider either as Lambert or Mercator
            if lat.mean()>30 or lat.mean()<-30:
                subplot_parameters = {"projection": ccrs.LambertConformal(central_longitude=ds.attrs["central_longitude"],
                                    central_latitude=ds.attrs["central_latitude"])}
            else:
                subplot_parameters = {"projection": ccrs.Mercator( central_longitude=ds.attrs["central_longitude"], 
                        min_latitude=lat.min(), max_latitude=lat.max() )}
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
            # cbar = plt.colorbar(temp_plot, ax=ax, orientation="vertical", label=myvar.attrs["units"])
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
    lat = ds["lat"]
    lon = ds["lon"]
    projection = ds.attrs["projection"]
    if "lambert" in projection.lower():
        subplot_parameters = {"projection": ccrs.LambertConformal(central_longitude=ds.attrs["central_longitude"], 
                    central_latitude=ds.attrs["central_latitude"], standard_parallels=ds.attrs["standard_parallels"],)}
    elif "mercator" in projection.lower():
        subplot_parameters = {"projection": ccrs.Mercator( central_longitude=ds.attrs["central_longitude"], 
                    min_latitude=lat.min(), max_latitude=lat.max() )}
    elif "polar" in projection.lower():
        if lat.mean() > 0:
            subplot_parameters = {"projection": ccrs.NorthPolarStereo( central_longitude=ds.attrs["central_longitude"] )}
        else:
            subplot_parameters = {"projection": ccrs.SouthPolarStereo( central_longitude=ds.attrs["central_longitude"] )}
    else: #consider either as Lambert or Mercator
        if lat.mean()>30 or lat.mean()<-30:
            subplot_parameters = {"projection": ccrs.LambertConformal(central_longitude=ds.attrs["central_longitude"],
                                central_latitude=ds.attrs["central_latitude"])}
        else:
            subplot_parameters = {"projection": ccrs.Mercator( central_longitude=ds.attrs["central_longitude"], 
                    min_latitude=lat.min(), max_latitude=lat.max() )}

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
    # cbar = plt.colorbar(temp_plot, ax=ax, orientation="vertical", label=myvar.attrs["units"])
    plt.colorbar(temp_plot, ax=ax, orientation="vertical", label=globals()[band].attrs["units"])
    plt.title("RGB Image from Brightness Temperatures")
    plt.axis("off")
    plt.savefig(postprocessingDir + "/" + "brightness_temperature_rgb.png")

if __name__ == "__main__":
    if not postprocessingEnabled:
        make_inputdata()
    else:
        postprocessing_directory_suffix = namelist["postprocessing"]['postprocessing_directory_suffix']
        postprocessingDir = wrfFileName+"_"+postprocessing_directory_suffix
        image_plot_enabled = namelist["postprocessing"]['image_plot_all_bands']
        rgb_plot_enabled = namelist["postprocessing"]['RGB_plot_brightness_temperature']['enabled']
        print("Postprocessing ...")
        #if not os.path.isdir(postprocessingDir):
        if not os.path.isdir(outputDirPath):
            print(f"Warning: The postprocessing directory ({postprocessingDir}) and RTTOV outputs are not availabe.")
            print("Disable postprocessing and run to make the profile files.")
            print("Exiting ..")
            exit()
        print(f"Converting the RTTOV output within the ({postprocessingDir}) directory to NetCDF files ..")
        make_netcdf()
        if image_plot_enabled:
            if os.path.isdir(postprocessingDir):
                print(f"Looking for the RTTOV NetCDF outputs in {postprocessingDir} ..")
                plot_png()
        if rgb_plot_enabled:
            if os.path.isdir(postprocessingDir):
                print(f"Looking for the RTTOV NetCDF outputs for brightness temperature in {postprocessingDir} ..")
                plot_rgb()
