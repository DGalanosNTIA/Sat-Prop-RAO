from sat_prop.core import wgs84ToEcef
from numpy import isclose

# input test case
e = 200  # meters
n = 100  # meters
u = 50   # meters
lat = 10 # degrees
lon = 12 # degrees
el = 55  # meters

# expected output
expected_x = 6144652
expected_y = 1306086
expected_z = 1100258

# function output
x, y, z = wgs84ToEcef(lat, lon, el)

# test
assert isclose(x, expected_x, rel_tol=1)
assert isclose(y, expected_y, rel_tol=1)
assert isclose(z, expected_z, rel_tol=1)

