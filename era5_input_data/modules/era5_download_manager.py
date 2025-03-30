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
    filename_level = "era5data_pressure_levels_"+str(year)+"_"+str(month)+"_"+str(day)+"_"+areaName+".nc"
    filename_surface = "era5data_surface_level_"+str(year)+"_"+str(month)+"_"+str(day)+"_"+areaName+".nc"
    c = cdsapi.Client()
    myData =     {
            'product_type': 'reanalysis',
            'variable': [
                'specific_humidity', 'temperature', 
            ],
            'pressure_level': [
                '1','2','3','5','7','10','20','30','50','70','100','125','150','175','200','225','250','300','350','400',
                '450','500','550','600','650','700','750','775','800','825','850','875','900','925','950','975','1000',
            ],
            'year': year,
            'month': month,
            'day': day,
            'time': hour,
            'area': [ north, west, south, east ],
            'format': 'netcdf',
        }
    print("  Downloading ERA5 data on the atmospheric pressure levels ..")
    print("-------------------------------------------------------------")
    #c.retrieve('reanalysis-era5-pressure-levels', myData, filename_level)
    c.retrieve('reanalysis-era5-pressure-levels', myData).download(filename_level)
    myData =     {
        'product_type': 'reanalysis',
        'variable': [
            '2m_temperature', '2m_dewpoint_temperature', 'surface_pressure', '10m_u_component_of_wind', '10m_v_component_of_wind',
            'skin_temperature', 'land_sea_mask', 'lake_cover', 'soil_type', 'geopotential', 'cloud_base_height', 'total_cloud_cover', 
        ],
        'year': year,
        'month': month,
        'day': day,
        'time': hour,
        'area': [ north, west, south, east ],
        'format': 'netcdf',
    }
    print("  Downloading ERA5 surface data ..")
    print("----------------------------------")
    #c.retrieve('reanalysis-era5-single-levels', myData, filename_surface)
    c.retrieve('reanalysis-era5-single-levels', myData).download(filename_surface)

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
