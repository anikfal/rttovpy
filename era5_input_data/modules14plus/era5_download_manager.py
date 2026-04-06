def main_dm():

    import importlib
    try:
       importlib.import_module("cdsapi")
    except:
       print("Warning: The Python module cdsapi is not installed.")
       print("Install it and run again.")
       print("Before running again, don't forget to set the CDS API key")
       print("Exiting ..")
       exit()

    import cdsapi
    import yaml
    with open('namelist_era5.yaml', 'r') as yaml_file:
        input_data = yaml.safe_load(yaml_file)
    year = input_data["time_of_simulation"]["year"]
    month = input_data["time_of_simulation"]["month"]
    day = input_data["time_of_simulation"]["day"]
    hour = input_data["time_of_simulation"]["hour"]
    check_integer(year)
    check_integer(month)
    check_integer(day)
    check_integer(hour)
    north = input_data["area_of_simulation"]["north_latitude"]
    south = input_data["area_of_simulation"]["south_latitude"]
    west  = input_data["area_of_simulation"]["west_longitude"]
    east = input_data["area_of_simulation"]["east_longitude"]
    areaName = input_data["area_of_simulation"]["domain_name"]
    check_float(north)
    check_float(south)
    check_float(west)
    check_float(east)
    filename_level = "era5data_pressure_levels_"+str(year)+"_"+str(month)+"_"+str(day)+"_"+areaName+".grib"
    filename_surface = "era5data_surface_level_"+str(year)+"_"+str(month)+"_"+str(day)+"_"+areaName+".grib"
    c = cdsapi.Client()
    
    #pressure level data
    myData =     {
            "class": "ea",
            "param": "130/133",
            "stream": "oper",
            "expver": "1",
            "levelist": "/".join(str(i) for i in range(1, 138)), #levels 1 to 137
            "levtype": "ml",
            "date": str(year)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2),
            'time': hour,
            'area': [ north, west, south, east ],
            "grid": [0.25, 0.25],
            "type": "an",
        }
    print("  Downloading ERA5 data on the atmospheric pressure levels ..")
    c.retrieve('reanalysis-era5-complete', myData).download(filename_level)
    print("-------------------------------------------------------------")
    
    #single level data
    myData =     {
        "class": "ea",
        "param": "134.128/167.128/168.128/165.128/166.128/235.128/172.128/228.128/43.128/129.128/164.128/228014.128",
        "stream": "oper",
        "expver": "1",
        "levtype": "sfc",
        "date": str(year)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2),
        'time': hour,
        'area': [ north, west, south, east ],
        "grid": [0.25, 0.25],
        "type": "an",
    }
    print("  Downloading ERA5 surface data ..")
    c.retrieve('reanalysis-era5-complete', myData).download(filename_surface)
    print("----------------------------------")

def check_integer(string_as_integer):
    try:
        int(string_as_integer)
    except ValueError:
        print("Warning: value", string_as_integer,  "in <input.yaml> is not an integer. Please check it out and correct it.")
        print("Exiting ..")
        exit()

def check_float(string_as_float):
    try:
        int(string_as_float)
    except ValueError:
        print("Warning: value", string_as_float,  "in <input.yaml> is not numeric. Please check it out and correct it.")
        print("Exiting ..")
        exit()

if __name__ == "__main__":
    main_dm()
