import era5_download_manager
import netCDF4 as nc
from glob import glob
import numpy as np
import math

# era5_download_manager.main_dm()

def specific_humidity_to_ppmv(specific_humidity):
    molar_mass_air = 28.97  # g/mol
    molar_mass_water_vapor = 18.02  # g/mol
    mixing_ratio = specific_humidity / (1 - specific_humidity)
    ppmv000 = mixing_ratio * (molar_mass_air / molar_mass_water_vapor) * 1e6
    ppmv = np.round(ppmv000).astype(int)
    return ppmv

def eSat(temperature):
    a1 = 611.21 #Pa
    a3_liquid = 17.502
    a4_liquid = 32.19 #K
    a3_ice = 22.587
    a4_ice = -0.7 #K
    T0 = 273.16
    if temperature>T0:
        ratio = (temperature - T0) / (temperature - a4_liquid)
        e_saturation = a1 * math.exp(a3_liquid * ratio)
        return e_saturation
    else:
        ratio = (temperature - T0) / (temperature - a4_ice)
        e_saturation = a1 * math.exp(a3_ice * ratio)
        return e_saturation
    
def qSat(temperature, pressure):
    epsilon = 0.622 #R_air/R_vapor
    saturation_vapor = (epsilon * eSat(temperature)) / (pressure - (1-epsilon)*eSat(temperature))
    return saturation_vapor


print(qSat(270, 50000))
print(specific_humidity_to_ppmv(qSat(270, 50000)))
exit()



ncFile_surface = nc.Dataset(glob("mydata_surface_*.nc")[0])
t2m = ncFile_surface["t2m"]
d2m = ncFile_surface["d2m"]
t2m = ncFile_surface["t2m"]
t2m = ncFile_surface["t2m"]



exit()

ncFile_level = nc.Dataset(glob("mydata_levels_*.nc")[0])
temperature000 = ncFile_level["t"]
qv = ncFile_level["q"]

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
        "1\n",
        "!\n",
    ]
    file_writer.writelines(header_lines)
file_append = open(profile_file, "a")
varShape = qv.shape
#for tt in range(varShape[0]):
for ii in range(varShape[2]): #latitudes
    for jj in range(varShape[3]): #longitude
        print(ii,jj)
        file_append.write("! --- Profile 1 ---\n!\n")
        for level in range(varShape[1]):
            pointValue = qv[0,level,ii,jj]
            file_append.write(str(pointValue)+'\n')
        exit()

file_append.close()

# with open(profile_file, "a") as file_appender:
#     file_appender.writelines(header_lines)

# print(qv[:5, 22, 22, 22])
