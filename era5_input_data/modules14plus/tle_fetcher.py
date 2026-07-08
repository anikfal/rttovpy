import os

# TLEs downloaded successfully are cached here, so a temporary server outage
# (e.g. CelesTrak returning HTTP 500) does not stop the run: the most recent
# cached TLE is used instead, with a warning about its age.
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tle_cache")

def _orbital_from_text(sat_name, tle_text):
    import tempfile
    from pyorbital.orbital import Orbital

    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp.write(tle_text)
    tmp.close()

    try:
        orb = Orbital(sat_name, tle_file=tmp.name)
        return orb
    finally:
        os.remove(tmp.name)

def _looks_like_tle(text):
    return any(line.startswith("1 ") for line in text.splitlines())

def _save_cache(cache_name, tle_text):
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(os.path.join(CACHE_DIR, cache_name + ".tle"), "w") as f:
        f.write(tle_text)

def _load_cache(cache_name):
    import time

    path = os.path.join(CACHE_DIR, cache_name + ".tle")
    if not os.path.isfile(path):
        return None, None
    age_days = (time.time() - os.path.getmtime(path)) / 86400.0
    with open(path) as f:
        return f.read(), age_days

def get_tle_celestrak(sat_name, sat_url):
    import requests

    try:
        response = requests.get(sat_url, timeout=10)
        response.raise_for_status()
        tle_text = response.text
        if not _looks_like_tle(tle_text):
            raise RuntimeError(
                f"CelesTrak returned no TLE data: {tle_text.strip()[:100]!r}"
            )
    except Exception as err:
        tle_text, age_days = _load_cache(sat_name)
        if tle_text is None:
            raise RuntimeError(
                f"CelesTrak request failed for {sat_name} and no cached TLE exists yet: {err}"
            )
        print(f"WARNING: CelesTrak request failed ({err}).")
        print(f"         Using cached TLE for {sat_name}, {age_days:.1f} days old.")
        return _orbital_from_text(sat_name, tle_text)

    _save_cache(sat_name, tle_text)
    return _orbital_from_text(sat_name, tle_text)

def get_tle_spacetrack_history(sat_name, norad_id, observation_time, username, password):
    import requests
    from datetime import timedelta

    login_url = "https://www.space-track.org/ajaxauth/login"

    start = observation_time - timedelta(days=1)
    end   = observation_time + timedelta(days=1)

    epoch = f"{start:%Y-%m-%d}--{end:%Y-%m-%d}"

    data_url = (
        "https://www.space-track.org/basicspacedata/query/"
        f"class/gp_history/NORAD_CAT_ID/{norad_id}/"
        f"EPOCH/{epoch}/orderby/EPOCH%20asc/format/tle"
    )

    # Historical TLEs are tied to the observation date, so cache per date.
    cache_name = f"{sat_name}_{observation_time:%Y%m%d}"

    try:
        session = requests.Session()
        login_response = session.post(login_url, data={
            "identity": username,
            "password": password
        }, timeout=10)
        if login_response.status_code != 200 or "Failed" in login_response.text:
            raise RuntimeError(
                "Space-Track login failed. Check space-track.org_username and "
                "space-track.org_password in the namelist."
            )

        response = session.get(data_url, timeout=10)
        response.raise_for_status()
        tle_text = response.text
        if not _looks_like_tle(tle_text):
            raise RuntimeError(
                f"Space-Track returned no TLE data for NORAD {norad_id}, epoch {epoch}"
            )
    except Exception as err:
        tle_text, age_days = _load_cache(cache_name)
        if tle_text is None:
            raise RuntimeError(
                f"Space-Track request failed for {sat_name} and no cached TLE exists yet: {err}"
            )
        print(f"WARNING: Space-Track request failed ({err}).")
        print(f"         Using cached historical TLE for {sat_name}, downloaded {age_days:.1f} days ago.")
        return _orbital_from_text(sat_name, tle_text)

    _save_cache(cache_name, tle_text)
    return _orbital_from_text(sat_name, tle_text)
