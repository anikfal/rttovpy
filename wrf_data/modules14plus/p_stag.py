import numpy as np
from wrf import getvar

def compute_p_stag(wrffile, observationIndex):
    # --- read variables ---
    P    = np.array(getvar(wrffile, "P",    timeidx=observationIndex))
    PB   = np.array(getvar(wrffile, "PB",   timeidx=observationIndex))
    PH   = np.array(getvar(wrffile, "PH",   timeidx=observationIndex))
    PHB  = np.array(getvar(wrffile, "PHB",  timeidx=observationIndex))
    PSFC = np.array(getvar(wrffile, "PSFC", timeidx=observationIndex))

    # --- full pressure ---
    P_total = P + PB

    # --- heights ---
    g = 9.81
    Z_stag = (PH + PHB) / g
    Z_mass = 0.5 * (Z_stag[:-1] + Z_stag[1:])

    nz, ny, nx = P_total.shape
    P_stag = np.zeros((nz + 1, ny, nx))

    # --- interpolation ---
    for i in range(ny):
        for j in range(nx):
            z_m = Z_mass[:, i, j]
            z_s = Z_stag[:, i, j]
            p_m = P_total[:, i, j]

            if z_m[0] > z_m[-1]:
                z_m = z_m[::-1]
                p_m = p_m[::-1]

            if z_s[0] > z_s[-1]:
                z_s = z_s[::-1]
                reversed_output = True
            else:
                reversed_output = False

            p_interp = np.exp(
                np.interp(z_s, z_m, np.log(p_m))
            )

            if reversed_output:
                p_interp = p_interp[::-1]

            P_stag[:, i, j] = p_interp

    # --- enforce ordering: TOA → surface ---
    P_stag = P_stag[::-1, :, :]

    # --- top extrapolation ---
    logp1 = np.log(P_stag[1])
    logp2 = np.log(P_stag[2])

    dz1 = Z_stag[-2] - Z_stag[-3]
    dz2 = Z_stag[-1] - Z_stag[-2]

    logp_top = logp1 + (logp1 - logp2) * (dz2 / dz1)
    P_stag[0] = np.exp(logp_top)

    # --- surface ---
    P_stag[-1] = PSFC

    return P_stag / 100.0  # hPa