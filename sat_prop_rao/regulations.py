'''
This module specifies limits on transmit power based on ITU recommendations.
'''

# TODO: Add in bandwidth of the Tx signal


def ITU_R_RA1631(frequency, antenna_type="dish"):
    '''
    :param frequency: central frequency of the RAO receiver in Hz
    :param antenna_type: antenna_type of antenna or threshold to determine limit offset or value (Options: dish, array, damage)
    :return: rx_epfd_limit
    '''
    # Limits are from ITU-R RA.769-2 and offsets/damage thresholds are from NRAO RFI Memo#153
    spectral_lines = [327, 1420, 1612, 1665, 4830, 14488,
                      22200, 23700, 43000, 48000, 88600, 150000, 220000, 265000]
    spectral_lines__Hz = [x * 1e06 for x in spectral_lines]
    assumed_bw = [10, 20, 20, 20, 50, 150, 250,
                  250, 500, 500, 1000, 1000, 1000, 1000]
    assumed_bw__Hz = [x * 1e03 for x in assumed_bw]
    pfd_limit__dBWperm2 = [-204, -196, -194, -194, -183, -
                           169, -162, -161, -153, -152, -148, -144, -139, -137]
    if antenna_type == "dish" or antenna_type == "array":
        # If the frequency falls above or below the max or min spectral lines, respectively, use the offset equation
        # based on the pfd_limit corresponding to the max or min spectral lines.
        if frequency < spectral_lines__Hz[0] - assumed_bw__Hz[0]/2:
            rx_epfd_limit = pfd_limit__dBWperm2[0] - \
                10/49 * frequency/1e9 + 46.2
            if antenna_type == "array":
                rx_epfd_limit = rx_epfd_limit + 10
            return rx_epfd_limit
        elif frequency > spectral_lines__Hz[len(spectral_lines__Hz)-1] + assumed_bw__Hz[len(spectral_lines__Hz)-1]/2:
            rx_epfd_limit = pfd_limit__dBWperm2[len(
                spectral_lines__Hz)-1] - 10 / 49 * frequency/1e9 + 46.2
            if antenna_type == "array":
                rx_epfd_limit = rx_epfd_limit + 10
            return rx_epfd_limit
        # Determine if the transmission falls on a spectral line.
        for i in range(len(spectral_lines__Hz)-1):
            # If the transmission falls on a spectral line, the stricter RA.769 Table 2 limit applies.
            # If not, there is still concern about LNA saturation.
            # NRAO RFI Memo 153 describes the required adjustment to RA.769 Table 2 threshold levels
            if spectral_lines__Hz[i] - assumed_bw__Hz[i] / 2 <= frequency <= spectral_lines__Hz[i] + assumed_bw__Hz[i] / 2:
                rx_epfd_limit = pfd_limit__dBWperm2[i]
                if antenna_type == "array":
                    rx_epfd_limit = rx_epfd_limit + 10
                return rx_epfd_limit
            elif i > 0 and \
                    spectral_lines__Hz[i-1] + assumed_bw__Hz[i-1]/2 < frequency < spectral_lines__Hz[i] - assumed_bw__Hz[i]/2:
                rx_epfd_limit = pfd_limit__dBWperm2[i -
                                                    1] - 10 / 49 * frequency/1e9 + 46.2
                if antenna_type == "array":
                    rx_epfd_limit = rx_epfd_limit + 10
                return rx_epfd_limit
    elif antenna_type == "damage":
        return -79
