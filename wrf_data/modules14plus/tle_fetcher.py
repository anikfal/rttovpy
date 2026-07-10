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

def _fetch_celestrak(sat_url):
    import requests

    response = requests.get(sat_url, timeout=10)
    response.raise_for_status()
    tle_text = response.text
    if not _looks_like_tle(tle_text):
        raise RuntimeError(
            f"CelesTrak returned no TLE data: {tle_text.strip()[:100]!r}"
        )
    return tle_text

def _fetch_tle_api(sat_name, norad_id):
    import requests

    response = requests.get(
        f"https://tle.ivanstanojevic.me/api/tle/{norad_id}", timeout=10,
        # The API closes the connection without responding when it sees the
        # default python-requests User-Agent, so an identifying one is sent.
        headers={"User-Agent": "rttovpy (satellite radiative transfer tool)"},
    )
    response.raise_for_status()
    data = response.json()
    # The API may report a different satellite name than the namelist uses
    # (e.g. "NOAA 20 (JPSS-1)"), so the title line is rewritten with sat_name
    # to keep pyorbital's name lookup working.
    tle_text = f"{sat_name}\n{data['line1']}\n{data['line2']}\n"
    if not _looks_like_tle(tle_text):
        raise RuntimeError(f"TLE API returned no TLE data for NORAD {norad_id}")
    return tle_text

def get_tle_celestrak(sat_name, sat_url, norad_id=None):
    sources = [("CelesTrak", lambda: _fetch_celestrak(sat_url))]
    if norad_id is not None:
        sources.append(
            ("TLE API (tle.ivanstanojevic.me)",
             lambda: _fetch_tle_api(sat_name, norad_id))
        )

    errors = []
    for source_name, fetch in sources:
        try:
            tle_text = fetch()
        except Exception as err:
            errors.append(f"{source_name}: {err}")
            continue
        if errors:
            print(f"WARNING: {'; '.join(errors)}")
            print(f"         TLE for {sat_name} downloaded from {source_name} instead.")
        _save_cache(sat_name, tle_text)
        return _orbital_from_text(sat_name, tle_text)

    all_errors = "; ".join(errors)
    tle_text, age_days = _load_cache(sat_name)
    if tle_text is None:
        raise RuntimeError(
            f"All TLE sources failed for {sat_name} and no cached TLE exists yet: {all_errors}"
        )
    print(f"WARNING: All TLE sources failed ({all_errors}).")
    print(f"         Using cached TLE for {sat_name}, {age_days:.1f} days old.")
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
