## Code for reading a WRF model output,and making RTTOV profile data based on the output data,
## as well as generating the shellscript for running the RTTOV model using the profile data.
## Author: Amirhossein Nikfal <https://github.com/anikfal>
###############################################################################
import argparse
def parse_args():
    parser = argparse.ArgumentParser(
        description="Run RTTOVpy or inspect coefficient wavelengths"
    )
    parser.add_argument(
        "--wavelength",
        type=str,
        help="Path to RTTOV coefficient file"
    )
    parser.add_argument(
        "--wrftime",
        type=str,
        help="Path to WRF output file to inspect its time range"
    )
    return parser.parse_args()
args = parse_args()
if args.wrftime:
    from modules14plus.rttov_utils import print_wrf_time_range
    print_wrf_time_range(args.wrftime)
    exit(0)
if args.wavelength:
    from modules14plus.rttov_utils import parse_rttov_wavenumbers, classify_channels
    wn_map = parse_rttov_wavenumbers(args.wavelength)
    info = classify_channels(wn_map)
    print("\nChannel | Wavenumber (cm⁻¹) | Wavelength (µm) | Solar")
    print("-" * 60)
    for ch in sorted(info):
        d = info[ch]
        print(f"{ch:7d} | {d['wavenumber']:16.3f} | {d['wavelength_um']:14.3f} | {d['requires_solar']}")
    exit(0)

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
with open('namelist_wrf.yaml', 'r') as yaml_file:
    namelist = yaml.safe_load(yaml_file)
rttovVersion = namelist["rttov_version"]

import sys
if __name__ == "__main__":
    try:
        if (rttovVersion == 12 or rttovVersion == 13):
            # executes on import → must be inside try
            from modules import run_rttov12
        else:
            from modules14plus import run_rttov14
            run_rttov14.main()
    except ValueError as e:
        print(f"\nError:\n{e}")
        sys.exit(1)