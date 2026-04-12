import numpy as np
import xarray as xr


def compute_era5_half_level_pressure(level_ds: xr.Dataset,
                                     sp: xr.DataArray) -> xr.DataArray:
    """
    Compute ERA5 L137 hybrid half-level pressure.

    Parameters
    ----------
    level_ds : xr.Dataset
        ERA5 model-level dataset opened with:
        backend_kwargs={"read_keys": ["pv", "NV"]}

    sp : xr.DataArray
        Surface pressure [Pa], dims (latitude, longitude)

    Returns
    -------
    xr.DataArray
        Half-level pressure [Pa]
        dims = ("half_level", "latitude", "longitude")
    """

    # ---- 1. Extract hybrid coefficients ----
    if "GRIB_pv" not in level_ds.t.attrs:
        raise RuntimeError(
            "Hybrid coefficients (pv) not found. "
            "Open GRIB with backend_kwargs={'read_keys':['pv','NV']}."
        )

    pv = np.asarray(level_ds.t.attrs["GRIB_pv"])

    n_half = pv.size // 2
    n_full = level_ds.sizes["hybrid"]   # <-- use sizes (no warning)

    if n_half != n_full + 1:
        raise RuntimeError(
            f"Inconsistent hybrid structure: "
            f"{n_full} full levels but {n_half} half levels."
        )

    a_half = pv[:n_half]
    b_half = pv[n_half:]

    # ---- 2. Compute half-level pressure ----
    # Convert sp to numpy array if needed
    sp_np = sp.values

    p_half = (
        a_half[:, None, None]
        + b_half[:, None, None] * sp_np
    )

    # ---- 3. Build xarray DataArray ----
    pressure_half = xr.DataArray(
        p_half,
        dims=("half_level", "latitude", "longitude"),
        coords={
            "half_level": np.arange(1, n_half + 1),
            "latitude": level_ds.latitude,
            "longitude": level_ds.longitude,
        },
        name="pressure_half",
        attrs={
            "units": "Pa",
            "long_name": "ERA5 hybrid half-level pressure"
        }
    )

    return pressure_half