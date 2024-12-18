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

TX_EFFICIENCY = .8
RX_EFFICIENCY = 0.72


def pointToPoint(frequency=3e9,
                 eirp=0,
                 txLatLonEl=[40, -150, 400e3],
                 rxLatLonEl=[45, -150, 0],
                 # pointing dir: east, north, up (ENU) unit vector
                 txDirENU=[0, 0, -1],
                 rxDirENU=[0, 0, 1],  # rx points up by default, tx points down

                 txAntPattern=antenna_pattern.rotationallySymmetric13dB,
                 rxAntPattern=antenna_pattern.rotationallySymmetric13dB,
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

    txAntGain = antenna_pattern.efficiencyAdjustment(txAntPattern(txTheta, 0), TX_EFFICIENCY)
    rxAntGain = antenna_pattern.efficiencyAdjustment(rxAntPattern(rxTheta, 0), RX_EFFICIENCY)

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

def testTxPFD(txLatLonEl = [0, 0, 400e3], tx_f_start = 1.990e9, tx_f_stop = 1.995e9, tx_target = [19.79, 0, 0],
              txAntPattern = antenna_pattern.ITU_S1528, propModel = propagation_model.freespace, txPowerIn__dBW = 20):
    radiusEarth__m = 6371e3
    alt__m = txLatLonEl[2]
    # map from lat/lon/el to earth-centered earth-fixed (ecef) coordinate system
    txPos = np.array(utl.wgs84ToEcef(*txLatLonEl))
    txTargetPos = np.array(utl.wgs84ToEcef(*tx_target))

    txToTxTargetVector = utl.unitVector(txTargetPos - txPos)

    # map from local heading ENU to ecef unit vector corresponding to heading
    txSubsatelliteVector = np.array(utl.enuToEcefHeading(*[0, 0, -1], *txLatLonEl))
    txPointingAngle__rad = utl.relativeAngle(txToTxTargetVector, txSubsatelliteVector)
    print("tx is pointing " + str(np.rad2deg(txPointingAngle__rad)) + " degrees off nadir.")
    satToTxTargetToEarthAngle__rad = pi-np.arcsin((radiusEarth__m+alt__m)/radiusEarth__m * np.sin(txPointingAngle__rad))
    print("The Earth center to satellite target to satellite angle is " + str(np.rad2deg(satToTxTargetToEarthAngle__rad)) + " degrees")
    satElAboveHorizon__rad = satToTxTargetToEarthAngle__rad-pi/2
    print("The satellite is " + str(np.rad2deg(satElAboveHorizon__rad)) + " degrees above the horizon.")
    print("The sum of the calculated internal angles is: " + str(np.rad2deg(satToTxTargetToEarthAngle__rad)+np.rad2deg(txPointingAngle__rad) + tx_target[0]-txLatLonEl[0]))
    if satElAboveHorizon__rad<=0:
        print("Satellite below horizon.")
        return None

    print("antenna 100% efficiency gain = " + str(txAntPattern(0,0)))
    txAntGain__dBi = txAntPattern(0, 0)
    # txAntGain__dBi = antenna_pattern.efficiencyAdjustment(txAntPattern(0, 0), TX_EFFICIENCY)
    print("antenna gain = " + str(txAntGain__dBi))
    # txAntCompensationGain = satGainAdjustment(txPointingAngle, txLatLonEl[2] / 1e3)


    min_distance__m = txLatLonEl[2]
    txChannelBW__MHz = 10 ** ((58 + 2.33) / 10) / 1e6  # ~1.08 MHz
    print("Tx Channel BW = " + str(txChannelBW__MHz))
    max_pfd = txPowerIn__dBW + txAntGain__dBi - 10*np.log10(4*pi*min_distance__m**2) - 10*np.log10(txChannelBW__MHz)
    target_pfd = -80
    available_padding__dB = max_pfd-target_pfd
    print("available padding: " + str(available_padding__dB))

    # calculate pathloss
    tx_f_c = (tx_f_start+tx_f_stop)/2
    distance__m = norm(txPos - txTargetPos)
    print("distance = " + str (distance__m))
    PL = propModel(distance__m, tx_f_c)
    print("Prop Loss = " + str(PL))

    txAntCompensationGain__dB = 10*np.log10(distance__m**2 / txLatLonEl[2]**2)
    print("initial Compensation Gain = " + str(txAntCompensationGain__dB))
    if txAntCompensationGain__dB > available_padding__dB:
        txAntCompensationGain__dB = available_padding__dB
    print("Final Compensation Gain = " + str(txAntCompensationGain__dB))
    print("Compensation Gain = " + str(txAntCompensationGain__dB))
    # power flux density


    pfd = txPowerIn__dBW + txAntGain__dBi - available_padding__dB + txAntCompensationGain__dB - 10*np.log10(4*pi*distance__m**2) - 10*np.log10(txChannelBW__MHz)

    print("Degrees above Horizon = " + str(np.rad2deg(satElAboveHorizon__rad)) + " pfd = " + str(pfd))

    return pfd

def satGainAdjustment(theta__rad, altitude__km):
    """
    To do the cdf calculation, we need to realistic adjusted gain of the satellite. The satellites change their gain
    to achieve their target epfd at the Earth's surface.  This adjustment needs to be considered in order to get an
    accurate epfd.
    From "A Framework for Assessing the Interference from NGSO Satellite Systems to a Radio Astronomy System"
    The compensation gain is dependent on the distance from the satellite to the targeted terminal and the altitude of
    the satellite.  In this case, we will consider the targeted terminal to be at the center of the boresight pointing
    direction so the gain compensation only depends on where it's pointing.
    """
    radius_of_Earth__km = 6371
    orbital_radius__km = radius_of_Earth__km + altitude__km
    """General Equation:  
    orbital_radius__km * np.cos(theta__rad)
    +/- np.sqrt(radius_of_Earth__km**2 - orbital_radius__km**2 * (np.sin(theta__rad))**2)"""
    d_sqrt_component__km = np.sqrt(radius_of_Earth__km**2 - orbital_radius__km**2 * (np.sin(theta__rad))**2)
    d_cos_component__km = orbital_radius__km * np.cos(theta__rad)
    if d_cos_component__km - d_sqrt_component__km < 0:
        d_terminal__km = d_cos_component__km + d_sqrt_component__km
    else:
        d_terminal__km = d_cos_component__km - d_sqrt_component__km

    compensatedGain = d_terminal__km**2 / altitude__km**2
    return compensatedGain

def testExclusionZone(frequency=3e9,
                      eirp=0,
                      txAntPattern=antenna_pattern.rotationallySymmetric13dB,
                      propModel=propagation_model.freespace
                      ):
    '''Developmental, unverified'''
    # EPFD limit defined by ITU-R RA1631 (todo: list additional relevant regs)
    epfd_limit = regulations.ITU_R_RA1631(frequency, antenna_type='dish')

    # N*N pixels in the heat map
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

    # build surface (include antennas here)
    X = np.zeros((N, N))
    for i in range(len(rx_pfd)):
        for j in range(len(txAntGain_arr)):
            X[i, j] = rx_pfd[i] + txAntGain_arr[j]

    # find contour line for EPFD limit in case0
    i_contour = np.zeros(N, dtype=int)
    for i in range(X.shape[0]):
        x_line = X[:, i]
        i_contour[i] = int(np.abs(x_line - epfd_limit).argmin())
    contour_case0 = y_axis[i_contour]
    
    # find contour line for EPFD limit in case0
    i_contour = np.zeros(N, dtype=int)
    for i in range(X.shape[0]):
        x_line = X[:, i]
        i_contour[i] = int(np.abs(x_line - (epfd_limit + 60)).argmin())
    contour_case1 = y_axis[i_contour]
    
    # find contour line for EPFD limit in case0
    i_contour = np.zeros(N, dtype=int)
    for i in range(X.shape[0]):
        x_line = X[:, i]
        i_contour[i] = int(np.abs(x_line - (epfd_limit + 100)).argmin())
    contour_case2 = y_axis[i_contour]

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
    plt.plot(x_axis, contour_case0, 'r')
    plt.plot(x_axis, contour_case1, 'r')
    plt.plot(x_axis, contour_case2, 'r')
    plt.gca().set_yticks(ytick_positions, ytick_labels)
    plt.xlabel('Offset Angle (degrees)')
    plt.ylabel('Separation Distance (km)')
    plt.title('Exclusion Zone Test')
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 15
    cbar.set_label('EPFD dB(W/(m^2))', rotation=270)
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
    # testExclusionZone()
    testTxPFD()

if __name__ == "__main__":
    main()
