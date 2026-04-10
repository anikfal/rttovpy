## Code for downloading ERA5 data,and making RTTOV profile data based on the downloaded data,
## As well as generating the shellscript for running the RTTOV model using the profile data
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################
import importlib
required_modules = ["yaml"]
for module in required_modules:
    try:
        importlib.import_module(module)
    except:
        print("Warning: The Python module", module, "is not installed.")
        print("Install it and run again.")
        print("Exiting ..")
        exit()
import yaml
with open('namelist_era5.yaml', 'r') as yaml_file:
    namelist = yaml.safe_load(yaml_file)
rttovVersion = namelist["rttov_version"]
if (rttovVersion == 12 or rttovVersion == 13):
    from modules import run_rttov12
else:
    from modules14plus import run_rttov14