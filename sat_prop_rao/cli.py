import click
import json
import sat_prop_rao.core as core
import sat_prop_rao.propagation_model as propagation_model
import sat_prop_rao.antenna_pattern as antenna_pattern


def parseFloatList(s): return [float(x) for x in s.split(',')]


antennaPatterns = {'isotropic': antenna_pattern.isotropic,
                   'basicRotationallySymmetric': antenna_pattern.basicRotationallySymmetric}
pathlossModels = {'freespace': propagation_model.freespace}


@click.command()
@click.option("--freq", default=1e9, help="Frequency (hz) of the RF signal.  [default: 1e9]", show_default=False)
@click.option("--eirp", default=0.0, help="Equivalent isotropic radiated power (EIRP) in dBm.  [default: 0 dBm]", show_default=False)
@click.option("--txPos", default='40,-105,400e3', help="Latitude (deg), longitude (deg), elevation (meters) of the transmitter.  [default: 40,-105,400e3]", show_default=False)
@click.option("--rxPos", default='40,-104,2', help="Latitude (deg), longitude (deg), elevation (meters) of the receiver.", show_default=True)
@click.option("--txAntENU", default='0,0,-1', help="The Array Normal of the transmit antenna in local east, north, up coordinates.", show_default=True)
@click.option("--rxAntENU", default='0,0,1', help="The Array Normal of the receive antenna in local east, north, up coordinates.", show_default=True)
@click.option("--txAntPattern", type=click.Choice(antennaPatterns), default='isotropic', help='Select antenna pattern for the transmitter.')
@click.option("--rxAntPattern", type=click.Choice(antennaPatterns), default='isotropic', help='Select antenna pattern for the receiver.')
@click.option("--propModel", type=click.Choice(pathlossModels), default='freespace', help='Propagation model.')
def pointToPoint(freq, eirp, txpos, rxpos, txantenu, rxantenu, txantpattern, rxantpattern, propmodel):
    """Estimates the point to point power flux density (PFD) in dB(W/m^2) from a satellite to a ground-based receiver.
       Output includes received signal strength (rx), power flux density (pfd), equivalent power flux density (epfd), and the epfd limit: 
       {"rx dBm": -170.36,
        "pfd dB(W/m^2)":-200.23,
        "epfd dB(W/m^2)":-180.24,
        "rx_epfd_limit dB(W/m^2)":-220.0}"""

    rx_dBm, rx_pfd, rx_epfd, rx_epfd_limit =\
        core.pointToPoint(frequency=freq,
                          eirp=eirp,
                          txLatLonEl=parseFloatList(txpos),
                          rxLatLonEl=parseFloatList(rxpos),
                          txDirENU=parseFloatList(txantenu),
                          rxDirENU=parseFloatList(rxantenu),

                          txAntPattern=antennaPatterns[txantpattern],
                          rxAntPattern=antennaPatterns[rxantpattern],
                          propModel=pathlossModels[propmodel]
                          )

    out = {}
    out['rx dBm'] = rx_dBm
    out['pfd dB(W/m^2)'] = rx_pfd
    out['epfd dB(W/m^2)'] = rx_epfd  # I think these units are wrong
    out['rx_epfd_limit dB(W/m^2)'] = rx_epfd_limit
    click.echo(json.dumps(out, indent=2))


@click.command()
@click.option("--freq", default=1e9, help="Frequency (hz) of the RF signal.  [default: 1e9]", show_default=False)
@click.option("--eirp", default=0.0, help="Equivalent isotropic radiated power (EIRP) in dBm.  [default: 0 dBm]", show_default=False)
@click.option("--txAntPattern", type=click.Choice(antennaPatterns), default='basicRotationallySymmetric', help='Select antenna pattern for the transmitter.')
@click.option("--propModel", type=click.Choice(pathlossModels), default='freespace', help='Propagation model.')
def exclusionZone(freq, eirp, txantpattern, propmodel):
    """Estimates the satellite exclusion zone around a typical radio astronomy observatory. Output shows the off-angle, distance contour line
       for the following cases:
         1. Receiver main beam pointing at satellite (worst case)
         2. Receiver -60 dB sidelobe pointing at satellite
         3. Receiver null -100 dB pointing at satellite
       This is currently a test function and is still a work in progress.
"""

    core.testExclusionZone(frequency=freq,
                               eirp=eirp,
                               txAntPattern=antennaPatterns[txantpattern],
                               propModel=pathlossModels[propmodel]
                               )

    click.echo("Todo: output Geo-json")