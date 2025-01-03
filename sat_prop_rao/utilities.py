'''
This module contains various utility functions to support vector math, unit conversion.
'''
import warnings
import numpy as np
from numpy.linalg import norm
from bisect import bisect_left
import random
from skyfield.api import load, Topos
from skyfield.iokit import parse_tle_file
from skyfield.toposlib import wgs84
from scipy.optimize import fsolve
from scipy.spatial import ConvexHull, Delaunay


def wgs84SurfaceLineIntersection(t, x0,y0,z0,vx,vy,vz):
    # WGS84 ellipsoid constants
    a = 6378137.0  # Semi-major axis (meters)
    b = 6356752.3142  # Semi-minor axis (meters)
    return ((x0 + t*vx)**2)/(a**2) + ((y0 + t*vy)**2)/(a**2) + ((z0 + t*vz)**2)/(b**2) - 1

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
    
    k = 100
    # north -> positive latitude (for small n)
    lat2 = lat + n / k

    # east -> positive longitude (for small e)
    lon2 = lon + e / k

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

class FieldOfViewOnEarthSurface():
    def __init__(self,observerLatLonEl,boundaryPrecisionDegrees=1):
        
        # observer position in ECEF
        observerPosition = np.array(wgs84ToEcef(*observerLatLonEl))
        
        # generate lat/lon lattice spanning earth
        latLimit = 90
        lonLimit = 180
        lats = np.linspace(-latLimit,latLimit,int((2*latLimit+1)/boundaryPrecisionDegrees))
        lons = np.linspace(-lonLimit,lonLimit,int((2*lonLimit+1)/boundaryPrecisionDegrees))
        
        # compute which lat,lon pairs are in and out of FOV
        latLonInFOV = []
        latLonOutFOV = []
        for lat in lats:
            for lon in lons:
                surfaceLatLonEl = [lat, lon, 1]
                surfacePosition = np.array(wgs84ToEcef(*surfaceLatLonEl))
                
                # vector from surface point to observer position
                vector = unitVector(observerPosition - surfacePosition)
                
                # check if vector intersects WGS84 ellipsoid
                initalGuess = 0 # corresponds to start of vector (surfacePosition)
                with warnings.catch_warnings() as caught_warnings:
                    warnings.filterwarnings("ignore", message="The iteration is not making good progress")
                    solution = fsolve(wgs84SurfaceLineIntersection,initalGuess,args=(*surfacePosition,*vector))
                
                if wgs84SurfaceLineIntersection(solution,*surfacePosition,*vector) > 1e-10:
                    latLonInFOV.append([lat,lon])  # the solution did not converge -> no intersection -> in FOV
                elif solution < 0:
                    latLonInFOV.append([lat,lon])  # the solution converged behind the start point -> in FOV
                else:
                    latLonOutFOV.append([lat,lon]) # the solution converged in front of start point -> out FOV
        
        # keeping these around for debugging purposes, if this function gets called a lot this should be deleted
        self._latLonInFOV = np.array(latLonInFOV)
        self._latLonOutFOV = np.array(latLonOutFOV)
        
        # convex hull (used to determine perimeter
        hull = ConvexHull(self._latLonInFOV)
        self.perimeter = self._latLonInFOV[hull.vertices]
        
        # used to determine if a point is within perimeter
        self._triangulate = Delaunay(self.perimeter)

        self.latBoundary = np.array([np.min(self.perimeter[:,0]), np.max(self.perimeter[:,0])])
        self.lonBoundary = np.array([np.min(self.perimeter[:,1]), np.max(self.perimeter[:,1])])

    def pointInFOV(self,latLon):
        return self._triangulate.find_simplex(latLon) >= 0
        
    def uniformSample(self,N):
        out = np.zeros((N,2))
        K = 0
        filled = False
        while not filled:
            batch = np.zeros((N,2))
            batch[:,0] = np.random.uniform(low=self.latBoundary[0],high=self.latBoundary[1],size=N)
            batch[:,1] = np.random.uniform(low=self.lonBoundary[0],high=self.lonBoundary[1],size=N)
            
            # points must be within FOV
            valid = batch[self.pointInFOV(batch)]
            k = valid.shape[0]
            
            if K+k <= N:
                out[K:K+k,:] = valid
            else:
                out[K:N,:] = valid[0:N-K]
                filled = True
            # try again until N points within FOV are generated
            K += k
        return out
    
    def uniformSampleECEF(self,N):
        x = self.uniformSample(N)
        out = np.zeros((x.shape[0],3))
        for i in range(x.shape[0]):
            out[i,:] = np.array(wgs84ToEcef(*x[i],0))
        return out

def dBmToWatts(x_dBm):
    return 10 ** ((x_dBm - 30) / 10)


def wattsTodBm(x_watts):
    return 10 * np.log10(x_watts) + 30

def satelliteMovement(altitude__km):
    """
    input: altitude in km
    output: circular velocity in rad/s
    For the cdf calculation we need to input a period of time that the measurement is being done and calculate the epfd
    during that time for a set of frequencies.  During that time, the satellites will move around.  We will consider all
    satellites moving generally from west to east but with a slight tilt for any satellites not on the equator.
     Their linear velocity is dependent on orbital mechanics.  The below is an approximation of these orbital mechanics.
    """
    gravitational_const__Nm2perkg2 = 6.6743*10**-11
    mass_of_Earth__kg = 5.9722*10**24
    radius_of_Earth__m = 6371000
    orbital_radius__m = radius_of_Earth__m + altitude__km*1000
    linear_velocity__mpers = np.sqrt(gravitational_const__Nm2perkg2 * mass_of_Earth__kg/orbital_radius__m)
    angular_velocity__radpers = 1/orbital_radius__m * linear_velocity__mpers

    return angular_velocity__radpers


def take_closest(my_list, my_number):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(my_list, my_number)
    if pos == 0:
        return my_list[0]
    if pos == len(my_list):
        return my_list[-1]
    before = my_list[pos - 1]
    after = my_list[pos]
    if after - my_number < my_number - before:
        return after, pos
    else:
        return before, pos - 1

def fieldOfViewToCells(target_solid_angle):
    elevation_step = np.sqrt(target_solid_angle)
    # The possible number of elevation rings that allows for even division in degrees
    el_step_values = [1, 2, 3, 5, 6, 9, 10, 15, 18, 30, 45, 90]
    # get the nearest step size from the eligible list
    elevation_step, el_index = take_closest(el_step_values, elevation_step)
    # The possible number of azimuth cells for even division in degrees
    # could also be found by a simple for loop with modulus, but since the eligible values are fixed this will be faster
    az_step_values = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360]
    min_el_angle = 0
    j = 0
    azimuthal_step, az_index = take_closest(az_step_values, elevation_step)
    az_list = []
    el_list = []
    sa_list = []
    while min_el_angle != 90:
        solid_angle_list = []
        for i in range(az_index, len(az_step_values)-1):
            solid_angle_list.append(np.sin(np.deg2rad(90 - (min_el_angle + elevation_step / 2)))
                                    * elevation_step * az_step_values[i])
        solid_angle = min(solid_angle_list, key=lambda x: abs(x - target_solid_angle))
        az_index = az_index + solid_angle_list.index(solid_angle)
        azimuthal_step = az_step_values[az_index]
        if min_el_angle == 0:
            g_mat = [[[0 for x in range(2)] for y in range(int(360 / azimuthal_step))] for z in range(int(90 / elevation_step))]
        for i in range(int(360 / azimuthal_step)):
            if min_el_angle != 87:
                # g_mat[j][i] = [min_el_angle, azimuthal_step * i, min_el_angle + elevation_step, azimuthal_step * (i+1)]
                g_mat[j][i] = [min_el_angle + elevation_step/2, azimuthal_step*(i+0.5)]
            else:
                # g_mat[j][i] = [min_el_angle, azimuthal_step * i, min_el_angle + elevation_step, 0]
                g_mat[j][i] = [min_el_angle + elevation_step/2, azimuthal_step*(i+0.5)]
        az_list.append(az_step_values[az_index])
        el_list.append(min_el_angle)
        sa_list.append(round(solid_angle, 2))
        j = j + 1
        min_el_angle = min_el_angle + elevation_step

    return az_list, el_list, sa_list, g_mat

def ecef_to_gps(x, y, z):
    """Convert ECEF coordinates back to GPS coordinates."""
    lat = np.degrees(np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2)))
    lon = np.degrees(np.arctan2(y, x))

    return lat, lon

def tle_import():
    max_days = 7.0 # download again once 7 days old
    name = 'starlink.txt' # custom filename to prevent overwrite for other data

    base = 'http://www.celestrak.org/NORAD/elements/gp.php'
    url = base + '?GROUP=starlink&FORMAT=tle'
    # TODO: Eventually we will want to change the below to save the date to the filename so we can go back and check old data.
    if not load.exists(name) or load.days_old(name) >= max_days:
        satellites = load.tle_file('http://www.celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle') # TLE data source
    else:
        ts = load.timescale()
        with load.open(name) as f:
            satellites = list(parse_tle_file(f, ts))

    # Set the initial time for the observation
    ts = load.timescale()
    t = ts.now() # set a specific (or random) time later
    min_lat = 0
    max_lat = 0
    for i in satellites:
        geocentric = i.at(t)
        sub_lat, sub_lon = wgs84.latlon_of(geocentric)

        if sub_lat.degrees < min_lat:
            min_lat = sub_lat.degrees
        elif sub_lat.degrees > max_lat:
            max_lat = sub_lat.degrees

    return min_lat, max_lat

def calculate_subsatellite_point(telescope_lat, telescope_lon, telescope_alt, azimuth, elevation, satellite_altitude):
    # Convert telescope coordinates to ECEF
    x0, y0, z0 = wgs84ToEcef(telescope_lat, telescope_lon, telescope_alt)

    # Convert azimuth and elevation to radians
    azimuth_rad = np.radians(azimuth)
    elevation_rad = np.radians(elevation)

    # Calculate the satellite position in ECEF coordinates
    slant_range = satellite_altitude  # Altitude of the satellite above the Earth
    x_s = x0 + slant_range * np.cos(elevation_rad) * np.cos(azimuth_rad)
    y_s = y0 + slant_range * np.cos(elevation_rad) * np.sin(azimuth_rad)
    z_s = z0 + slant_range * np.sin(elevation_rad)

    # Convert back to GPS coordinates
    subsatellite_lat, subsatellite_lon = ecef_to_gps(x_s, y_s, z_s)
    return subsatellite_lat, subsatellite_lon, satellite_altitude

def randomSatPosSelection(ras_lat, ras_lon, ras_alt, los_cell_matrix, num_satellites, sat_alt):
    min_lat, max_lat = tle_import()
    sat_el_list = []
    sat_az_list = []
    distance_list = []
    ras_loc = wgs84.latlon(ras_lat, ras_lon)
    for i in range(num_satellites):
        sat_ind_row = random.randint(0, int(len(los_cell_matrix[0]) / 4)-1)
        sat_ind_col = random.randint(0, int(len(los_cell_matrix[1]))-1)
        print("Max row = ", int(len(los_cell_matrix[0]) / 4)-1)
        print("Max col = ", int(len(los_cell_matrix[1]))-1)
        print("sat_ind_row = ", sat_ind_row)
        print("sat_ind_col = ", sat_ind_col)
        sat_el, sat_az = los_cell_matrix[sat_ind_row][sat_ind_col]
        sat_altitude = random.choice(sat_alt)
        sat_lat, sat_lon, sat_altitude = calculate_subsatellite_point(ras_lat, ras_lon, ras_alt, sat_az, sat_el, sat_altitude)
        while sat_lat > max_lat or sat_lat < min_lat:
            sat_ind_row = random.randint(0, int(len(los_cell_matrix[0]) / 4)-1)
            sat_ind_col = random.randint(0, int(len(los_cell_matrix[1]))-1)
            print("sat_ind_row2 = ", sat_ind_row)
            print("sat_ind_col2 = ", sat_ind_col)
            sat_el, sat_az = los_cell_matrix[sat_ind_row][sat_ind_col]
            sat_altitude = random.choice(sat_alt)
            sat_lat, sat_lon, sat_altitude = calculate_subsatellite_point(ras_lat, ras_lon, ras_alt, sat_az, sat_el, sat_altitude)

        sat_el_list.append(sat_el)
        sat_az_list.append(sat_az)

        ts = load.timescale()
        t = ts.now()  # set a specific (or random) time later
        satellite_pos = wgs84.latlon(sat_lat, sat_lon, sat_altitude)
        difference = satellite_pos - ras_loc
        topocentric = difference.at(t)
        alt, az, distance = topocentric.altaz()
        print("distance: ", distance.km)
        distance = float(distance.km)
        distance_list.append(distance)
    return sat_el_list, sat_az_list, distance_list

def unwantedEmissionPowerRASBand(f_c__Hz, B_n, atten_bw__kHz = 4, ):
    """
    From NTIA Manual of Regulations and Procedures for Federal Radio Frequency Management
    www.ntia.gov/sites/default/files/2023-03/complete_manual_january_2022_revision.pdf page 424
    NOTE: Not using this one at this time as the SCS satellites are only held to the FCC unwanted emissions limits
          described in the unwantedEmissionsLimitsFCC function.
    NOTE2: If this function is used, we should implement it similarly to the way the FCC one was done in the
           unwantedEmissionsLimitsFCC function
    f_c__Hz: center frequency of NGSO Transmitting channel
    B_n: necessary bandwidth of NGSO transmitting channel
    TODO: Annex J of the NTIA Manual of Regulations and Procedures for Federal Radio Frequency Management has the actual
          equations for B_n (likely use J.3.6 pg 782). We should implement those at some point. For now, we can just use
          the channel BW.
    atten_bw__kHz: This determines the roll-off based on the applicable BW of the channel.  Normally 4 kHz, but
                   increases to 1MHz for center frequencies above 15GHz
    """
    f_spurious = 10*B_n
    for i in range(0, f_spurious, f_spurious/1000):
        if i/B_n <= .5:
            psd__dB = 0
        elif i/B_n <=10:
            psd__dB = -40*np.log10(2*i/B_n) - 8
        else:
            psd__dB = -60

def unwantedEmissionLimitsFCC(RAS_res_BW__kHz, RAS_meas_f_c__MHz, RAS_meas_span__MHz,
                              scs_f_c__MHz, scs_channel_bw__MHz, scs_p_t__W):
    """
    From www.ecfr.gov/current/title-47/chapter-I/subchapter-B/part-25/subpart-C/section-25.202
    Inputs:
    RAS_res_BW__kHz:  radio astronomy system resolution bandwidth in kHz
    RAS_meas_f_c__MHz:  radio astronomy system measurement center frequency in MHz
    RAS_meas_span__MHz:  radio astronomy system measurement frequency span in MHz
    scs_f_c__MHz:  SCS satellite transmission center frequency in MHz
    scs_channel_bw__MHz:  SCS satellite transmission channel bandwidth in MHz
    scs_p_t__W:  SCS satellite transmission center frequency power in Watts (power leaving antenna)
    Output:
    tx_spectral_mask_dict__dB:  spectral mask dictionary describing the attenuation of the transmitted signal over the
                               frequency span of the RAS measurement with the frequency steps as the dictionary keys.
    TODO: Add in aggregate limit
    TODO: Create overall spectral mask based on bandwidth limitations of the RAS hardware instead of measurement setup
          as driving the receiver LNA into the non-linear region depends on the total power at the LNA, not just the
          power in the measurement setup.
    """
    atten_bw__kHz = 4
    tx_spectral_mask_dict__dB = {}
    for f__kHz in range(RAS_meas_f_c__MHz*10**3 - (RAS_meas_span__MHz*10**3)/2,
                        RAS_meas_f_c__MHz*10**3 + (RAS_meas_span__MHz*10**3)/2, RAS_res_BW__kHz):
        if f__kHz <= scs_f_c__MHz + 0.5*scs_channel_bw__MHz:
            tx_spectral_mask_dict__dB.update({str(f__kHz): 10*np.log10(scs_p_t__W) - 0})
        elif f__kHz <= scs_f_c__MHz + scs_channel_bw__MHz:
            tx_spectral_mask_dict__dB.update({str(f__kHz): 10*np.log10(scs_p_t__W) - 25})
        elif f__kHz <= scs_f_c__MHz + 2.5*scs_channel_bw__MHz:
            tx_spectral_mask_dict__dB.update({str(f__kHz): 10*np.log10(scs_p_t__W) - 35})
        else:
            tx_spectral_mask_dict__dB.update({str(f__kHz): 10*np.log10(scs_p_t__W) - (43 + 10*np.log10(scs_p_t__W))})

    return tx_spectral_mask_dict__dB

def receiverMask():
    """
    Placeholder for when/if we add in the receiver mask via code.  It may be best to do it solely through JSON in which
    case this code will pull the JSON file in and return the signal attenuation dictionary similar to
    unwantedEmissionsLimitsFCC.
    """

def epfdByFreq(rx_f_start, rx_f_stop):
    """
    Placeholder for a function to get the total EPFD over the receiving band including the out of band emissions from
    the transmitter and any filtering done by the receiver.
    """

if __name__ == '__main__':
    ras_lat = 38.4264
    ras_lon = -79.8370
    ras_alt = 807.43
    sat_alts = [340000, 550000, 1100000]

    az, el, sa, pos_mat = fieldOfViewToCells(9)
    print("Az = ", az)
    print("El = ", el)
    print("SA = ", sa)
    print("Position Matrix: ")
    for i in pos_mat:
        for j in i:
            print(j, end=' ')
        print()
    sat_ind_row = random.randint(0, int(len(pos_mat[0])/4))
    sat_ind_col = random.randint(0, int(len(pos_mat[1])))
    print("pos mat rows: ", int(len(pos_mat[0])/4))
    print("pos mat cols: ", int(len(pos_mat[1])))
    print("sat_col: ", sat_ind_col)
    print("sat_row: ", sat_ind_row)
    sat_el, sat_az = pos_mat[sat_ind_row][sat_ind_col]
    print("az: " + str(sat_az) + ", el: " + str(sat_el))
    sat_alt = random.choice(sat_alts)
    sat_lat, sat_lon, sat_alt = calculate_subsatellite_point(ras_lat, ras_lon, ras_alt, sat_az, sat_el, sat_alt)
    print("Sat Pos: " + str(sat_lat) + ", " + str(sat_lon) + " elev: " + str(sat_alt))
    min_lat, max_lat = tle_import()
    print("The maximum latitude Starlink is found at is: ", max_lat)
    print("The minimum latitude Starlink is found at is: ", min_lat)
    sat_el, sat_az, sat_distances = randomSatPosSelection(ras_lat, ras_lon, ras_alt, pos_mat, 30, sat_alts)
    print("Satellite Elevations: ", sat_el)
    print("Satellite Azimuths: ", sat_az)
    print("Satellite Distances: ", sat_distances)
