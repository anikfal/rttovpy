from datetime import datetime, timezone

def choose_tle_source(sim_time, high_precision):
    now = datetime.now(timezone.utc)
    delta = abs((now - sim_time).total_seconds())

    if not high_precision:
        return "celestrak"

    # threshold ~ 2 days (tunable)
    if delta < 2 * 86400:
        return "celestrak"
    else:
        return "spacetrack_history"