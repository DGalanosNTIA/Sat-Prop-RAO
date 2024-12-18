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
    gain0 = np.where(phi <= phi_m, g_max - 2.5*10**-3 *
                     (diameter/wavelength*phi)**2, 0)
    phi = np.where(phi <= phi_m, 200, phi)
    gain1 = np.where(phi <= phi_r, g_one, 0)
    phi = np.where(phi <= phi_r, 200, phi)
    gain2 = np.where(phi <= 10, 29 - 25 * np.log10(phi), 0)
    phi = np.where(phi <= 10, 200, phi)
    gain3 = np.where(phi <= 34.1, 34-30*np.log10(phi), 0)
    phi = np.where(phi <= 34.1, 200, phi)
    gain4 = np.where(phi <= 80, -12, 0)
    phi = np.where(phi <= 80, 200, phi)
    gain5 = np.where(phi <= 120, -7, 0)
    phi = np.where(phi <= 120, 200, phi)
    gain6 = np.where(phi <= 180, -12, 0)

    gain = gain0+gain1+gain2+gain3+gain4+gain5+gain6
    return gain

def ITU_S1528(az__rad, el__rad, max_gain__dBi = 38, z=1, L_N__dB = -60, L_F__dB = -38):
    # Assuming: Rotationally symmetric and determined piece-wise by the angle psi from boresight
    psi__degrees = np.rad2deg(np.sqrt((az__rad ** 2 + el__rad ** 2) / 2))
    # z = 1 for circular beams
    # one-half the 3 dB beamwidth in the plane of interest
    satPsi_b__degrees = (np.sqrt(4 * np.pi / (10 ** (max_gain__dBi / 10))) * 180 / np.pi) / 2
    a = 2.58*np.sqrt(1-.06*np.log10(z))
    b = 6.32
    alpha = 1.5
    Y = b * satPsi_b__degrees * 10**(0.04 * (max_gain__dBi + L_N__dB + L_F__dB))
    X = max_gain__dBi + L_N__dB + 25 * np.log10(z)
    L_B__dB_tmp = 15+L_N__dB+.25*max_gain__dBi+5*np.log10(z)
    if L_B__dB_tmp > 0:
        L_B__dB = L_B__dB_tmp
    else:
        L_B__dB = 0

    # where phi is less than ### replace with the equation in that angle range otherwise set to 0
    gain0__dBi = np.where(psi__degrees <= a * satPsi_b__degrees, max_gain__dBi - 3 *
                     (psi__degrees/satPsi_b__degrees)**alpha, 0)
    psi__degrees = np.where(psi__degrees <= a * satPsi_b__degrees, 200, psi__degrees)
    gain1__dBi = np.where(psi__degrees <= 0.5 * b * satPsi_b__degrees, max_gain__dBi + L_N__dB + 20*np.log10(z), 0)
    psi__degrees = np.where(psi__degrees <= 0.5 * b * satPsi_b__degrees, 200, psi__degrees)
    gain2__dBi = np.where(psi__degrees <= b * satPsi_b__degrees, max_gain__dBi + L_N__dB, 0)
    psi__degrees = np.where(psi__degrees <= b * satPsi_b__degrees, 200, psi__degrees)
    gain3__dBi = np.where(psi__degrees <= Y, X - 25 * np.log10(psi__degrees), 0)
    psi__degrees = np.where(psi__degrees <= Y, 200, psi__degrees)
    gain4__dBi = np.where(psi__degrees <= 90, L_F__dB, 0)
    psi__degrees = np.where(psi__degrees <= 90, 200, psi__degrees)
    gain5__dBi = np.where(psi__degrees <= 180, L_B__dB, 0)

    gain__dBi = gain0__dBi+gain1__dBi+gain2__dBi+gain3__dBi+gain4__dBi+gain5__dBi
    return gain__dBi

def efficiencyAdjustment(gain, efficiency):
    trueGain = gain * efficiency
    return trueGain



# utilities
#def rotationallySymmetricInterpolation(az, el, theta, gaindBi):
    #return 1


#def loadJson(az, el, filename):
    #return 1
