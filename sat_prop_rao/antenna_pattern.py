'''
This module contains several common antenna patterns as identified by the ITU
as well as utilities to import antenna pattern files. When including an
application-specific antenna pattern, consider using one of the utilities or
construct your own using the following format.

antennaPatternFunction:
Returns Antenna gain in dBi 

        Parameters:
                az (int): Azimuth angle in radians
                el (int): Elevation angle in radians
                freq (float): parameter to support freq dependent antenna patterns in Hz
                *args: Additional positional arguments. These arguments are stored in a tuple (specify as parameter if necessary).
                **kwargs: Additional keyword arguments. These arguments are stored in a dictionary (specify as parameter if necessary).

        Returns:
                antGain (float): dBi antenna gain
'''
import numpy as np
from scipy import constants

# Antenna Pattern examples:
def isotropic(az, el, freq=0):
    return 0*az*el

def rotationallySymmetric13dB(az, el, freq=0): return \
    20*np.log10(np.abs(np.sinc(np.sqrt(az**2 + el**2)/2)))

def ITU_RA1631(az, el, freq=1e9, diameter=10*(constants.c/1e9)):
    # Rotationally symmetric and determined piece-wise by the angle phi from boresight
    phi = np.rad2deg(np.sqrt((az**2 + el**2)/2))
    wavelength = constants.c / freq
    g_max = 20*np.log10(diameter/wavelength) + 20*np.log10(np.pi)
    g_one = -1 + 15*np.log10(diameter/wavelength)
    phi_m = 20*wavelength/diameter * np.sqrt(g_max - g_one)
    phi_r = 15.85 * (diameter/wavelength)**-0.6

    # where phi is less than ### replace with the equation in that angle range otherwise set to 0
    gain0 = np.where(phi < phi_m, g_max - 2.5*10**-3 *
                     (diameter/wavelength*phi)**2, 0)
    phi = np.where(phi < phi_m, 200, phi)
    gain1 = np.where(phi < phi_r, g_one, 0)
    phi = np.where(phi < phi_r, 200, phi)
    gain2 = np.where(phi < 10, 29 - 25 * np.log10(phi), 0)
    phi = np.where(phi < 10, 200, phi)
    gain3 = np.where(phi < 34.1, 34-30*np.log10(phi), 0)
    phi = np.where(phi < 34.1, 200, phi)
    gain4 = np.where(phi < 80, -12, 0)
    phi = np.where(phi < 80, 200, phi)
    gain5 = np.where(phi < 120, -7, 0)
    phi = np.where(phi < 120, 200, phi)
    gain6 = np.where(phi <= 180, -12, 0)

    gain = gain0+gain1+gain2+gain3+gain4+gain5+gain6
    return gain

# utilities
#def rotationallySymmetricInterpolation(az, el, theta, gaindBi):
    #return 1


#def loadJson(az, el, filename):
    #return 1
