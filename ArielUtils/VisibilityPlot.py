import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from ArielUtils.Functions import *

def plot_sky_view(data, 
                  title,
                  plot_ecliptic=True, 
                  cmap="hot", 
                  save_path="Visibility Plots/", 
                  long_pts=360, 
                  dpi=300, 
                  size=1, 
                  marker='o', 
                  label='Stars', 
                  marker_color='blue', 
                  projection='aitoff'):
    lat_pts = long_pts//2
    
    # preparing data
    starRA = np.radians(data['Star RA'])
    starRA = starRA.apply(lambda x: x-2*pi if x>pi else x) # ensures that RA values are between -pi and pi.
    starDEC = np.radians(data['Star Dec'])

    # coordinates
    ra = np.linspace(-pi, pi, long_pts)
    dec = np.linspace(-pi/2, pi/2, lat_pts)
    
    # preparing ecliptic
    sun_path = pd.DataFrame(np.array([[80, 172, 264, 355], 
                                      [0, pi/2, pi, -pi/2], 
                                      [0, np.radians(23.5), 0, np.radians(-23.5)]]).transpose(), 
                            columns=['Day of year', 'Ra', 'Dec'])
    ecliptic = sun_path['Dec'][1]*np.sin(ra)
    if plot_ecliptic:
        eclip = "_withEcliptic"
    else:
        eclip = ""

    # preparing data for the pcolormesh to be plotted
    Lon, Lat = np.meshgrid(ra, dec)

    # density
    if long_pts == 360:
        rho = np.genfromtxt('ArielUtils/VisibilityDensity_360.csv', delimiter=',')
    elif long_pts == 720:
        rho = np.genfromtxt('ArielUtils/VisibilityDensity_720.csv', delimiter=',')
    else:
        rho = np.zeros(np.shape(Lon)) # array that will contain all the values for the pcolormesh (higher value = more visible throughout year)
        for n in range(len(ra)):
            for i, row in enumerate(dec):
                for j, col in enumerate(ra):
                    if dist_angle(row,col, ecliptic[n],ra[n]) > np.radians(MIN_ANGLE_TOWARDS_SUN) and dist_angle(row, col, ecliptic[n]+pi, ra[n]) > np.radians(MIN_ANGLE_AWAY_SUN):
                        rho[i][j] += 1
        rho /= long_pts # normalizes values such that they are all between 0 and 1. 

    # Preparing the plot
    plt.figure(figsize=(10,5))
    plt.subplot(projection = projection)
    plt.grid(True, alpha=.6)
    
    # plotting the pcolormesh and its color bar
    plt.pcolormesh(Lon, Lat, rho, cmap=cmap)
    plt.colorbar(label="Fraction of the year where\na given location is visible")

    # plotting ecliptic path
    if plot_ecliptic:
        plt.plot(ra, ecliptic, label='Ecliptic', color='blue', linewidth=1) # plots ecliptic
        
    # plotting data
    plt.scatter(starRA, starDEC, s=size, marker=marker, label=label, color=marker_color) # plots the stars of interest

    plt.xticks([-pi/2, 0, pi/2], ['18h', '0h', '6h'], color='white')
    plt.yticks([-pi/3, -pi/6, 0, pi/6, pi/3])

    # adding title and labels
    plt.title(title, pad=20)
    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.legend(loc='upper right')

    # saving figure to the "Visibility Plots" folder
    plt.savefig(f"{save_path}{title}_CMAP={cmap}{eclip}.png", bbox_inches='tight', dpi=dpi)
    plt.show()