# Function to convert to Julian Date from a datetime object
def datetime_to_jd(dt):
    # Julian Date formula for converting datetime to JD
    a = (14 - dt.month) // 12
    y = dt.year + 4800 - a
    m = dt.month + 12 * a - 3
    jd = dt.day + ((153 * m + 2) // 5) + 365 * y + (y // 4) - (y // 100) + (y // 400) - 32045
    jd += (dt.hour - 12) / 24.0 + dt.minute / 1440.0 + dt.second / 86400.0
    return jd
