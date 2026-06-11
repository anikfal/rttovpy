import numpy as np
from wrf import getvar

# --- read variables from WRF ---
P    = getvar(wrffile, "P",    timeidx=observationIndex)
PB   = getvar(wrffile, "PB",   timeidx=observationIndex)
PH   = getvar(wrffile, "PH",   timeidx=observationIndex)
PHB  = getvar(wrffile, "PHB",  timeidx=observationIndex)
PSFC = getvar(wrffile, "PSFC", timeidx=observationIndex)

# --- convert to numpy arrays (avoid xarray issues) ---
P    = np.array(P)
PB   = np.array(PB)
PH   = np.array(PH)
PHB  = np.array(PHB)
PSFC = np.array(PSFC)

# --- full-level pressure (Pa) ---
P_total = P + PB

# --- geopotential height ---
g = 9.81
Z_stag = (PH + PHB) / g
Z_mass = 0.5 * (Z_stag[:-1] + Z_stag[1:])

# --- allocate ---
nz, ny, nx = P_total.shape
P_stag = np.zeros((nz + 1, ny, nx))

# --- logarithmic interpolation ---
for i in range(ny):
    for j in range(nx):
        z_m = Z_mass[:, i, j]
        z_s = Z_stag[:, i, j]
        p_m = P_total[:, i, j]

        # ensure increasing order
        if z_m[0] > z_m[-1]:
            z_m = z_m[::-1]
            p_m = p_m[::-1]

        if z_s[0] > z_s[-1]:
            z_s = z_s[::-1]
            reversed_output = True
        else:
            reversed_output = False

        p_interp = np.exp(
            np.interp(
                z_s,
                z_m,
                np.log(p_m)
            )
        )

        if reversed_output:
            p_interp = p_interp[::-1]

        P_stag[:, i, j] = p_interp

# --- fix vertical ordering (top → surface) ---
P_stag = P_stag[::-1, :, :]

# --- Option B: improve top level (log extrapolation) ---
logp1 = np.log(P_stag[1, :, :])
logp2 = np.log(P_stag[2, :, :])

dz1 = Z_stag[-2, :, :] - Z_stag[-3, :, :]
dz2 = Z_stag[-1, :, :] - Z_stag[-2, :, :]

logp_top = logp1 + (logp1 - logp2) * (dz2 / dz1)
P_stag[0, :, :] = np.exp(logp_top)

# --- surface (must enforce) ---
P_stag[-1, :, :] = PSFC

# --- convert to hPa ---
P_stag = P_stag / 100.0