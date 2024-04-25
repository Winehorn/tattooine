from astropy.coordinates import SkyCoord
from sunpy.coordinates import get_body_heliographic_stonyhurst
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np

REAL_RADIUS = True

obstime = Time('2022-02-26T19:10:00')
planet_list = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
planet_coord = [get_body_heliographic_stonyhurst(this_planet, time=obstime) for this_planet in planet_list]

fig = plt.figure()
ax1 = plt.subplot(1, 1, 1, projection='polar')

ax1.yaxis.grid(False)

def get_real_radius(coord):
    # Use ln to scale outer planet's radii.
    if this_coord.radius.value > 2:
        return np.log(this_coord.radius.value)
    else:
        return this_coord.radius.value

for this_planet, this_coord in zip(planet_list, planet_coord):
    if REAL_RADIUS:
        radius = get_real_radius(this_coord)
    else:
        radius = planet_list.index(this_planet) + 1
    plt.plot(np.deg2rad(this_coord.lon), radius, 'o', label=this_planet)
    
    # Add orbit.
    orbit = plt.Circle((0.0, 0.0), radius, transform=ax1.transData._b, fill=False)
    ax1.add_artist(orbit)
    print(f'Planet: {this_planet} | Radius: {this_coord.radius.value} | Displayed Radius: {radius}')

if REAL_RADIUS:
    inner_outer_border = plt.Circle((0.0, 0.0), planet_coord[3].radius.value + (np.log(planet_coord[4].radius.value) - planet_coord[3].radius.value) / 2, transform=ax1.transData._b, fill=False, linestyle='--')
    ax1.add_artist(inner_outer_border)
#plt.legend(loc='lower left')
plt.savefig('tattoo.pdf')


