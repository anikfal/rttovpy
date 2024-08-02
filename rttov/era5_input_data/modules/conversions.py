import numpy as np

def surface_humidity(t2m_dew, p2m):
    a1 = 611.21 #Pa
    a3_liquid = 17.502
    a4_liquid = 32.19 #K
    a3_ice = 22.587
    a4_ice = -0.7 #K
    T0 = 273.16 #K
    epsilon = 0.622 #R_air/R_vapor
    ratio = (t2m_dew - T0) / (t2m_dew - a4_liquid)
    e_vapor_liquid = a1 * np.exp(a3_liquid * ratio)
    ratio = (t2m_dew - T0) / (t2m_dew - a4_ice)
    e_vapor_ice = a1 * np.exp(a3_ice * ratio)
    mask_zero_liquid = np.where(t2m_dew>T0, 1, 0)
    mask_zero_ice = np.where(mask_zero_liquid==1, 0, 1)
    e_vapor_liquid_mask = e_vapor_liquid * mask_zero_liquid
    e_vapor_ice_mask = e_vapor_ice * mask_zero_ice
    e_vapor = e_vapor_liquid_mask + e_vapor_ice_mask
    specific_humidity = (epsilon * e_vapor) / (p2m - (1-epsilon) * e_vapor)
    return specific_humidity

def short2float(ncVar):
    ncVar = ncVar.scale_factor * ncVar + ncVar.add_offset
    return ncVar

def specific_humidity_to_ppmv(specific_humidity):
    epsilon = 0.622 #R_air/R_vapor
    mixing_ratio = specific_humidity / (1 - specific_humidity)
    ppmv000 = mixing_ratio * epsilon * 1e6
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
        e_saturation = a1 * np.exp(a3_liquid * ratio)
        return e_saturation
    else:
        ratio = (temperature - T0) / (temperature - a4_ice)
        e_saturation = a1 * np.exp(a3_ice * ratio)
        return e_saturation
    
def qSat(temperature, pressure):
    epsilon = 0.622 #R_air/R_vapor
    saturation_vapor = (epsilon * eSat(temperature)) / (pressure - (1-epsilon)*eSat(temperature))
    return saturation_vapor