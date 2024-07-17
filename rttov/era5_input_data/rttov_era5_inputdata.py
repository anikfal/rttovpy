import era5_download_manager
import netCDF4 as nc
from glob import glob
import numpy as np

# era5_download_manager.main_dm()

ncfile_level = nc.Dataset(glob("mydata_levels_*.nc")[0])
temperature000 = ncfile_level["t"]
ncfile_level = nc.Dataset(glob("mydata_levels_*.nc")[0])
qv = ncfile_level["q"]


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
varshape = qv.shape
#for tt in range(varshape[0]):
for ii in range(varshape[2]): #latitudes
    for jj in range(varshape[3]): #longitude
        print(ii,jj)
        file_append.write("! --- Profile 1 ---\n!\n")
        for kk in range(varshape[1]):
            file_append.write(str(qv[0,kk,ii,jj])+'\n')
        exit()

file_append.close()

# with open(profile_file, "a") as file_appender:
#     file_appender.writelines(header_lines)

# print(qv[:5, 22, 22, 22])
