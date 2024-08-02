def main_dm():
    import cdsapi
    import yaml
    with open('input.yaml', 'r') as yaml_file:
        input_data = yaml.safe_load(yaml_file)
    years = read_check_comma(input_data, "time", "years")
    months = read_check_comma(input_data, "time", "months")
    days = read_check_comma(input_data, "time", "days")
    hours000 = read_check_comma(input_data, "time", "hours")
    for element in years:
        check_integer(element)
    for element in months:
        check_integer(element)
    for element in days:
        check_integer(element)
    for element in hours000:
        check_integer(element)
    north = input_data["area"]["north_latitude"]
    south = input_data["area"]["south_latitude"]
    west  = input_data["area"]["west_longitude"]
    east = input_data["area"]["east_longitude"]
    areaName = input_data["area"]["name"]
    check_float(north)
    check_float(south)
    check_float(west)
    check_float(east)
    my_year = years[0]
    my_month = months[0]
    my_day = days[0]
    filename_level = "era5data_pressure_levels_"+my_year+"_"+my_month+"_"+my_day+"_"+areaName+".nc"
    filename_surface = "era5data_surface_level_"+my_year+"_"+my_month+"_"+my_day+"_"+areaName+".nc"
    hours = [hour.zfill(2)+":00" for hour in hours000]
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
            'year': years,
            'month': months,
            'day': days,
            'time': hours,
            'area': [ north, west, south, east ],
            'format': 'netcdf',
        }
    print("  Downloading ERA5 data on the atmospheric pressure levels ..")
    print("-------------------------------------------------------------")
    # c.retrieve('reanalysis-era5-pressure-levels', myData, filename_level)
    myData =     {
        'product_type': 'reanalysis',
        'variable': [
            '2m_temperature', '2m_dewpoint_temperature', 'surface_pressure', '10m_u_component_of_wind', '10m_v_component_of_wind',
            'skin_temperature', 'land_sea_mask', 'lake_cover', 'soil_type', 'geopotential', 'cloud_base_height', 'total_cloud_cover', 
        ],
        'year': years,
        'month': months,
        'day': days,
        'time': hours,
        'area': [ north, west, south, east ],
        'format': 'netcdf',
    }
    print("  Downloading ERA5 surface data ..")
    print("----------------------------------")
    c.retrieve('reanalysis-era5-single-levels', myData, filename_surface)

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

def read_check_comma(inputfile, param, var):
    value = inputfile[param][var]
    if isinstance(value, int):
        return [str(value)]
    else:
        value = ''.join(value.split())
        if value.endswith(","):
            value = value[:-1]
        return value.split(",")

if __name__ == "__main__":
    main_dm()