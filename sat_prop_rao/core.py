'''
This module contains the primary processing code to compute point-to-point
and exclusion zone estimates for satellite to ground station RF propagation.
These are available as command line tools:
sat-prop-ptp
sat-prop-exclusion
'''
import sat_prop_rao.antenna_pattern as antenna_pattern
import sat_prop_rao.propagation_model as propagation_model
import sat_prop_rao.regulations as regulations
import sat_prop_rao.utilities as utl
import numpy as np
from numpy.linalg import norm
from scipy.constants import pi
from scipy.constants import speed_of_light as c
from scipy.spatial.transform import Rotation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt


def pointToPoint(frequency=3e9,
                 eirp=0,
                 txLatLonEl=[40, -150, 400e3],
                 rxLatLonEl=[45, -150, 0],
                 # pointing dir: east, north, up (ENU) unit vector
                 txDirENU=[0, 0, -1],
                 rxDirENU=[0, 0, 1],  # rx points up by default, tx points down

                 txAntPattern=antenna_pattern.basicRotationallySymmetric,
                 rxAntPattern=antenna_pattern.basicRotationallySymmetric,
                 propModel=propagation_model.freespace,
                 ):
    # need to check that tx/rx are within a relevant distance of one another (e.g. remove case: opposite sides of the world)

    # map from lat/lon/el to earth-centered earth-fixed (ecef) coordinate system
    txPos = np.array(utl.wgs84ToEcef(*txLatLonEl))
    rxPos = np.array(utl.wgs84ToEcef(*rxLatLonEl))

    # map from local heading ENU to ecef unit vector corresponding to heading
    txHeading = np.array(utl.enuToEcefHeading(*txDirENU, *txLatLonEl))
    rxHeading = np.array(utl.enuToEcefHeading(*rxDirENU, *rxLatLonEl))

    # extract beam pattern gain (currently only for az/el symmetric antenna)
    txToRxHeading = utl.unitVector(rxPos - txPos)
    txTheta = utl.relativeAngle(txToRxHeading, txHeading)
    rxTheta = utl.relativeAngle(-txToRxHeading, rxHeading)

    txAntGain = txAntPattern(txTheta, 0)
    rxAntGain = rxAntPattern(rxTheta, 0)

    # calculate pathloss
    distance = norm(txPos - rxPos)
    PL = propModel(distance, frequency)

    # power received
    rx_dBm = eirp - PL + txAntGain + rxAntGain

    # power flux density
    rx_pfd = rx_dBm - rxAntGain + 20 * \
        np.log10(frequency) - 20*np.log10(c) + 10*np.log10(4*pi)

    # equivalent power flux density
    rx_epfd = rx_pfd + rxAntGain

    # EPFD limit defined by ITU-R RA1631
    rx_epfd_limit = regulations.ITU_R_RA1631(frequency)

    return rx_dBm, rx_pfd, rx_epfd, rx_epfd_limit


def testExclusionZone(frequency=3e9,
                      eirp=0,
                      txAntPattern=antenna_pattern.basicRotationallySymmetric,
                      propModel=propagation_model.freespace
                      ):
    '''Developmental, unverified'''
    # EPFD limit defined by ITU-R RA1631 (todo: list additional relevant regs)
    epfd_limit = regulations.ITU_R_RA1631(frequency, antenna_type='dish')

    # span separation distance
    N = 1000

    # min/max distance in logspace to map (todo: make arguments)
    a = 10*np.log10(1e3)
    b = 10*np.log10(100e6)

    # plot axes values
    x_axis = np.linspace(-180, 180, N)
    y_axis = np.linspace(0, 1, N)

    # array of distances and corresponding pathloss values
    r_arr = 10**(np.linspace(a, b, N)/10)
    PL_arr = propModel(r_arr, frequency)

    # power received (w/o antennas considered)
    rx_dBm = eirp - PL_arr

    # power flux density (w/o antennas considered)
    rx_pfd = rx_dBm + 20 * \
        np.log10(frequency) - 20*np.log10(c) + 10*np.log10(4*pi)

    # span offset angle
    theta_arr = np.linspace(-pi, pi, N)
    txAntGain_arr = txAntPattern(theta_arr, 0)

    # build surface
    X = np.zeros((N, N))
    for i in range(len(rx_pfd)):
        for j in range(len(txAntGain_arr)):
            X[i, j] = rx_pfd[i] + txAntGain_arr[j]

    # find contour line for EPFD limit
    i_contour = np.zeros(N, dtype=int)
    for i in range(X.shape[0]):
        x_line = X[:, i]
        i_contour[i] = int(np.abs(x_line - epfd_limit).argmin())
    contour = y_axis[i_contour]

    # imshow doesn't seem to handle very large numbers on the axes well in this case
    # so I have to map everything to [0,1] and back out where to place the tick marks
    # need to verify this
    ytick_labels = np.array(
        [10, 100, 1000, 10000, 100000, 100000, 10000000], dtype=int)
    ytick_positions = []
    for i in range(len(ytick_labels)):
        j = np.abs(ytick_labels[i]*1e3 - r_arr).argmin()
        ytick_positions.append(y_axis[j])

    plt.figure(figsize=(10, 6))
    plt.imshow(X, extent=(-180, 180, 0, 1), aspect='auto')
    plt.plot(x_axis, contour, 'r')
    plt.gca().set_yticks(ytick_positions, ytick_labels)
    plt.xlabel('Offset Angle (degrees)')
    plt.ylabel('Separation Distance (km)')
    plt.title('Exclusion Zone Test')
    plt.colorbar()
    plt.show()

    # todo:
    # Considering three cases for RAO antenna gain:
    # 1. Receiver main beam pointing at satellite (worst case)
    # 2. Receiver -60 dB sidelobe pointing at satellite
    # 3. Receiver null -100 dB pointing at satellite
    # rxAntGainOptions = np.array([0, -60, -100])

    # for k in range(len(rxAntGainOptions)):
    #    rxAntGain = rxAntGainOptions[k]


def main():
    # print(pointToPoint())
    testExclusionZone()


if __name__ == "__main__":
    main()
