'''
This module contains various utility functions to support vector math, unit conversion.
'''
import numpy as np
from numpy.linalg import norm


def wgs84ToEcef(lat, lon, el):
    # WGS84 ellipsoid constants
    a = 6378137.0  # Semi-major axis (meters)
    b = 6356752.3142  # Semi-minor axis (meters)
    e_sq = 1 - (b ** 2 / a ** 2)  # First eccentricity squared

    # convert latitude and longitude from degrees to radians
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    # calculate N
    N = a / np.sqrt(1 - e_sq * np.sin(lat_rad) ** 2)

    # calculate ECEF coordinates
    X = (N + el) * np.cos(lat_rad) * np.cos(lon_rad)
    Y = (N + el) * np.cos(lat_rad) * np.sin(lon_rad)
    Z = ((b ** 2 / a ** 2) * N + el) * np.sin(lat_rad)

    return X, Y, Z


def enuToEcefHeading(e, n, u, lat, lon, el):
    # only valid for small ENU, used as a pointing vector

    # get ecef of reference
    x_ref, y_ref, z_ref = wgs84ToEcef(lat, lon, el)

    # force enu to be a unit vector
    e, n, u = unitVector(np.array([e, n, u]))

    # north -> positive latitude (for small n)
    lat2 = lat + n / 1e6

    # east -> positive longitude (for small e)
    lon2 = lon + e / 1e6

    # up and el add
    el2 = el + u

    # approximate solution
    x_rough, y_rough, z_rough = wgs84ToEcef(lat2, lon2, el2)

    # create a unit vector corresponding to ecef heading
    u = unitVector(
        np.array([x_rough - x_ref, y_rough - y_ref, z_rough - z_ref]))

    return u[0], u[1], u[2]


def unitVector(u):
    return u / norm(u)


def relativeAngle(u, v):
    # ensures angles are unitary
    u = unitVector(u)
    v = unitVector(v)
    return np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))


def dBmToWatts(x_dBm):
    return 10 ** ((x_dBm - 30) / 10)


def wattsTodBm(x_watts):
    return 10 * np.log10(x_watts) + 30
