
'''
This module contains relevant pathloss functions for sattelite to ground propagation.
Additional pathloss functions may be added in the form of:

propagationLossFunction:
Returns pathloss in dB 

        Parameters:
                distance (float): distance in meters
                frequency (float): frequency in Hz
                *args: Additional positional arguments. These arguments are stored in a tuple (specify as parameter if necessary).
                **kwargs: Additional keyword arguments. These arguments are stored in a dictionary (specify as parameter if necessary).

        Returns:
                pathloss (float): dB pathloss
'''
import numpy as np
from scipy.constants import pi
from scipy.constants import speed_of_light as c


def freespace(distance, frequency): return \
    20*np.log10(distance) + 20*np.log10(frequency) + 20*np.log10(4*pi/c)

# todo: include IF77
