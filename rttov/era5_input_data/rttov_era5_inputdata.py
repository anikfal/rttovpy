import era5_download_manager
import netCDF4 as nc
from glob import glob
import numpy as np

# era5_download_manager.main_dm()

ncfile_level = nc.Dataset(glob("mydata_levels_*.nc")[0])
temperature000 = ncfile_level["t"]
ncfile_level = nc.Dataset(glob("mydata_levels_*.nc")[0])
qv = ncfile_level["q"]


profile_filename = "myprofile.txt"
with open(profile_filename, "w") as file_writer:
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
        "2\n",
        "!\n",
    ]
    file_writer.writelines(header_lines)


# with open(profile_filename, "a") as file_appender:
#     file_appender.writelines(header_lines)

# print(qv[:5, 22, 22, 22])
