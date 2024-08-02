from modules import era5_download_manager
from modules.conversions import surface_humidity
import yaml
import netCDF4 as nc
from glob import glob
import numpy as np




with open('input.yaml', 'r') as yaml_file:
    inputFile = yaml.safe_load(yaml_file)
satIndex = inputFile["satellite_specifications"]["satellite_name_index"]
if satIndex > 76:
    print("Warning: Satellite name index must be between 1 to 76. Please look inside <satellite_names.yaml>.")
    print("exiting ..")
    exit()
angleEnable = inputFile["satellite_specifications"]["user_defined_angles"]['enable']

print("angleEnable:", angleEnable)

with open('satellite_names.yaml', 'r') as yaml_file:
    satNameFile = yaml.safe_load(yaml_file)
print("sat index:", satIndex, satNameFile[satIndex] )

if(angleEnable):
    print("yesssssssssssss")

exit()

# era5_download_manager.main_dm()


ncFile_surface = nc.Dataset(glob("era5data_surface_level_*.nc")[0])
t2m = ncFile_surface["t2m"]
d2m = ncFile_surface["d2m"]
p2m = ncFile_surface["sp"]
u10 = ncFile_surface["u10"]
v10 = ncFile_surface["v10"]
skinT = ncFile_surface["skt"]
landSeaMask = np.array(ncFile_surface["lsm"]).round()
lakeCover = ncFile_surface["cl"]
soilType = ncFile_surface["slt"]
geopotential = ncFile_surface["z"]
cloudBase = ncFile_surface["cbh"]
cloudCover = ncFile_surface["tcc"]
q2m = surface_humidity(t2m[:], p2m[:])
lat = ncFile_surface['latitude']
lon = ncFile_surface['longitude']


ncFile_level = nc.Dataset(glob("era5data_pressure_levels_*.nc")[0])
temperature000 = ncFile_level["t"]
qv = ncFile_level["q"]
tempLevel = ncFile_level["t"]
pressureLevels = ncFile_level['level']

profile_file = "myprofile.txt"
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
varShape = qv.shape
#for tt in range(varShape[0]):
tt = 0
for ii in range(varShape[2]): #latitudesTemperature profile (K)
    for jj in range(varShape[3]): #longitude
        print("ii:", ii, "  jj:", jj)
        file_append.write("! --- Profile 1 ---\n")
        subHead = [
                "!\n",
                "! Pressure levels (hPa)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        for level in range(varShape[1]):
            file_append.write(str(pressureLevels[level])+'\n')
        subHead = [
                "!\n",
                "! Temperature profile (K)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        for level in range(varShape[1]):
            pointValue = tempLevel[tt,level,ii,jj]
            file_append.write(str(pointValue)+'\n')
        subHead = [
                "!\n",
                "! Water vapour profile (kg/kg)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        for level in range(varShape[1]):
            pointValue = qv[0,level,ii,jj]
            file_append.write(str(pointValue)+'\n')
        subHead = [
                "!\n",
                "! Near-surface variables:\n"
                "!  2m T (K)    2m q (ppmv) 2m p (hPa) 10m wind u (m/s)  10m wind v (m/s)  wind fetch (m)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        near_surface_vars = [t2m[tt, ii, jj], q2m[tt, ii, jj], p2m[tt, ii, jj]/100,\
                              u10[tt, ii, jj], v10[tt, ii, jj], 100000]
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
        skinVars = [skinT[tt, ii, jj], 35.0, 3.0, 5.0, 15.0, 0.1, 0.3]
        skinVars_2line = ' '.join(map(str, skinVars))
        file_append.write(skinVars_2line)
        subHead = [
                "\n!\n",
                "! Surface type (0=land, 1=sea, 2=sea-ice) and water type (0=fresh, 1=ocean)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        surfaceType = [landSeaMask[tt, ii, jj], 1]
        surfaceType_2line = ' '.join(map(str, surfaceType))
        file_append.write(surfaceType_2line)
        subHead = [
                "\n!\n",
                " Elevation (km), latitude and longitude (degrees)\n"
                "!\n",
            ]
        for line in subHead:
            file_append.write(line)
        elevation = [geopotential[tt, ii, jj]/9.81, lat[jj], lon[ii]]
        elevation_2line = ' '.join(map(str, elevation))
        file_append.write(elevation_2line)
        exit()

file_append.close()

# with open(profile_file, "a") as file_appender:
#     file_appender.writelines(header_lines)

# print(qv[:5, 22, 22, 22])
