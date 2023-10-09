import numpy as np

# Sun constraints
MIN_ANGLE_TOWARDS_SUN = 70 # in degrees
MIN_ANGLE_AWAY_SUN = 60 # in degrees

# Telescope constraints
SLEWRATE = np.radians(4.5) # radians per minute
SETTLE_TIME = 5 # minutes

# Mission constrants
MISSION_DURATION = 0.5 *365*24*60 # minutes

# Calibration constraints
SHORT_CALIBRATION_DURATION = 1 *60 # minutes
SHORT_CALIBRATION_FREQUENCY = 36 *60 # minutes
LONG_CALIBRATION_DURATION = 6 *60 # minutes
LONG_CALIBRATION_FREQUENCY = 30 *24*60 # minutes

# Baseline parameters
BASELINE_RATIO = 0.25 # ratio of baseline-to-event-duration on either side of the event. 
BASELINE_DURATION = 30 # baseline duration, on either side of the event

# Original Reference time
T_0 = 0 # MJD calendar (November 17ty 1858)