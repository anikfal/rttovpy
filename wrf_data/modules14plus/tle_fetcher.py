def get_tle_celestrak(sat_name, sat_url):
    import requests
    import tempfile
    import os
    from pyorbital.orbital import Orbital

    tle_text = requests.get(sat_url, timeout=10).text

    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp.write(tle_text)
    tmp.close()

    try:
        orb = Orbital(sat_name, tle_file=tmp.name)
        return orb
    finally:
        os.remove(tmp.name)

def get_tle_spacetrack_history(sat_name, norad_id, observation_time, username, password):
    import requests
    import tempfile
    import os
    from datetime import timedelta
    from pyorbital.orbital import Orbital

    login_url = "https://www.space-track.org/ajaxauth/login"

    start = observation_time - timedelta(days=1)
    end   = observation_time + timedelta(days=1)

    epoch = f"{start:%Y-%m-%d}--{end:%Y-%m-%d}"

    data_url = (
        "https://www.space-track.org/basicspacedata/query/"
        f"class/gp_history/NORAD_CAT_ID/{norad_id}/"
        f"EPOCH/{epoch}/orderby/EPOCH%20asc/format/tle"
    )

    session = requests.Session()
    session.post(login_url, data={
        "identity": username,
        "password": password
    })

    tle_text = session.get(data_url, timeout=10).text

    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp.write(tle_text)
    tmp.close()

    try:
        orb = Orbital(sat_name, tle_file=tmp.name)
        return orb
    finally:
        os.remove(tmp.name)