# Functions for ARIEL scheduling and visualizations.

# # # importing useful libraries # # #
from numpy import arccos, sin, cos, pi
import pandas as pd

from astropy.time import Time
from astropy.coordinates import get_sun

from ArielUtils.Constants import *
# # # # # #


# # # FUNCTIONS # # #

# Distance between two points on the celestial sphere.
def dist_angle(dec1:float, 
               ra1:float, 
               dec2:float, 
               ra2:float)->float:
    """Returns the distance (in radians) between two points on the celestial sphere.

    Args:
        dec1 (float): Declination of first point (in radians).
        ra1 (float): Right ascention of first point (in radians).
        dec2 (float): Declination of second point (in radians).
        ra2 (float): Right ascention of second point (in radians).

    Returns:
        float: Distance between two given points on the celestial sphere, in radians.
    """
    ra_angle = cos(ra2 - ra1)
    return arccos( 0.5 * (cos(dec2-dec1)*(ra_angle+1) + cos(dec2 + dec1)*(ra_angle-1) ))


# Finds every visibe target at a given time
def find_visible_targets(targets:pd.DataFrame, 
                         time:Time=Time.now(), 
                         sun_angle:float=MIN_ANGLE_TOWARDS_SUN, 
                         opp_sun_angle:float=MIN_ANGLE_AWAY_SUN)->pd.DataFrame:
    """Returns a data frame of all the targets that are visible at a given time. 

    Args:
        targets (pd.DataFrame): Data frame containing the targets.
        time (Time, optional): given time, as Time object. Defaults to time.now().
        sun_angle (float, optional): Limiting angle towards Sun, in degrees. Defaults to MIN_ANGLE_TOWARDS_SUN.
        opp_sun_angle (float, optional): Limiting angle away from Sun, in degrees. Defaults to MIN_ANGLE_AWAY_SUN.

    Returns:
        pd.DataFrame: Data frame containing all the targets visible at a given time. 
    """
    # initializing the list of visible targets and the Sun's coordinates at the given time
    viewable_targets = []
    sunRA = np.radians(get_sun(time).ra.value) # Sun's right ascension, in radians
    sunDEC = np.radians(get_sun(time).dec.value) # Sun's declination, in radians
    
    # iterating through the targets data frame. _ just means that the index will not be used in the loop
    for _, target in targets.iterrows():
        distance = dist_angle(np.radians(target['Star Dec']), np.radians(target['Star RA']), sunDEC, sunRA)
        if distance >= np.radians(sun_angle) and distance <= pi - np.radians(opp_sun_angle): # checks if a target is visible at the given time
            viewable_targets.append(target) # adding target to list since it is visible at the given time
    return pd.DataFrame(viewable_targets)



# Finds the closest target to the current pointing position.
def find_closest_target(targets:pd.DataFrame, 
                        currentRA:float=0, 
                        currentDEC:float=0)->pd.DataFrame:
    """Returns a data frame containing the nearest target.

    Args:
        targets (pd.DataFrame): Data frame containing targets
        currentRA (float, optional): Current right ascension. Defaults to 0.
        currentDEC (float, optional): Current declination. Defaults to 0.

    Returns:
        pd.DataFrame: Data frame containing the target nearest current location
    """
    # initializing necessary variables
    closest_target = None 
    closest_distance = None
    
    # iterating through target data frame. _ just means that the index will not be used in the loop.
    for _, target in targets.iterrows():
        distance = dist_angle(np.radians(currentDEC), np.radians(currentRA), np.radians(target['Star Dec']), np.radians(target['Star RA']))
        if closest_target is None or distance < closest_distance: # checking if the target is closer than the previously closest target.
            closest_target = target
            closest_distance = distance
    return pd.DataFrame(closest_target).transpose()



# Finds the total slewtime, given a distance to cover.
def find_slewtime_minutes(dec1:float, 
                          ra1:float, 
                          dec2:float, 
                          ra2:float, 
                          slewrate:float=SLEWRATE)->float:
    """Returns the slewtime in minutes given the current position and the final/target position.

    Args:
        dec1 (float): initial declination, in radians.
        ra1 (float): initial right ascension, in radians.
        dec2 (float): target declination, in radians.
        ra2 (float): target right ascension, in radians.
        slewrate (float, optional): telescope slew rate. Defaults to SLEWRATE.

    Returns:
        float: slew time, in minutes.
    """
    distance = dist_angle(dec1, ra1, dec2, ra2)
    return (distance / slewrate)


# Checks if calibration must be done
def check_calibration(time_elapsed, frequency, error, time_of_previous_calibration):
    if (time_elapsed - time_of_previous_calibration) <= frequency-error:
        return False
    diff = time_elapsed % frequency
    if diff <= error or frequency-diff <= error:
        return True
    return False


# Fitness function
def fitness(distance:float, 
            event_duration:float, 
            orbital_period:float, 
            time_till_event:float, 
            quality_metric:float=1, 
            slewrate:float=SLEWRATE, 
            settle_time:float=SETTLE_TIME, 
            baseline_duration:float=BASELINE_RATIO):
    wait_time = time_till_event - (distance/slewrate + settle_time + baseline_duration/2)
    if wait_time < 0:
        F = 0 # setting fitness to 0 if wait time is negative
    else:
        F = orbital_period * quality_metric / (wait_time+1)
    return F, wait_time

# Sorts targets by fitness
def sort_by_fitness(targets:pd.DataFrame, 
                    time:Time=Time.now(),
                    currentRA:float=0, 
                    currentDEC:float=0)->list:
    targets_copy = targets.copy(deep=True)
    time.format = 'mjd'
    for index, target in targets_copy.iterrows():
        dist = dist_angle(currentDEC, currentRA, target['Star Dec'], target['Star RA'])
        period = target['Planet Period [days]']
        t_diff = time.value - target['T_ref']
        time_till_transit = (period - t_diff) % period
        fit_val, wait_time = fitness(dist, 
                          target['Transit Duration [s]']/60,
                          period*24*60, 
                          time_till_transit*24*60
                          )
        targets_copy.at[index, 'Fitness'] = fit_val
        targets_copy.at[index, 'Wait Time'] = wait_time
    targets_copy.sort_values(by=['Fitness'], ascending=False, inplace=True)
    return targets_copy
# # # # # #