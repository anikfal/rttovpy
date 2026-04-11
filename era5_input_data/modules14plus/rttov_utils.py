import re
import numpy as np

def parse_rttov_wavenumbers(coef_file):
    """
    Extract channel -> wavenumber mapping from RTTOV coef file.
    Returns dict: {channel_number: wavenumber_cm-1}
    """
    channel_wavenumbers = {}

    in_filter_block = False

    with open(coef_file, "r") as f:
        for line in f:
            line = line.strip()

            if "FILTER_FUNCTIONS" in line:
                in_filter_block = True
                continue

            if in_filter_block:
                if line.startswith("! ---") or line.startswith("FUNDAMENTAL_CONSTANTS"):
                    break

                if not line or line.startswith("!"):
                    continue

                parts = line.split()

                if len(parts) >= 3:
                    try:
                        ch = int(parts[0])
                        wavenumber = float(parts[2])
                        channel_wavenumbers[ch] = wavenumber
                    except:
                        continue

    return channel_wavenumbers

def classify_channels(channel_wavenumbers):
    """
    Returns dict:
    {
        ch: {
            "wavenumber": ...,
            "wavelength_um": ...,
            "requires_solar": True/False
        }
    }
    """
    result = {}

    for ch, wn in channel_wavenumbers.items():
        wavelength = 10000.0 / wn  # µm

        # robust rule
        requires_solar = wavelength < 4.0

        result[ch] = {
            "wavenumber": wn,
            "wavelength_um": wavelength,
            "requires_solar": requires_solar
        }

    return result

def requires_solar_for_channels(coef_file, selected_channels):
    channel_wavenumbers = parse_rttov_wavenumbers(coef_file)
    channel_info = classify_channels(channel_wavenumbers)

    solar_channels = []
    for ch in selected_channels:
        if ch not in channel_info:
            raise ValueError(f"Channel {ch} not found in coef file")

        if channel_info[ch]["requires_solar"]:
            solar_channels.append(ch)

    return len(solar_channels) > 0, solar_channels, channel_info