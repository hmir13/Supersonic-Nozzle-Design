from StandardAtmosphere import StandardAtmosphere
from Nozzle import Nozzle

import matplotlib.pyplot as plt
import numpy as np

# parameters [SI Units]
altitude = 10000    # [m]
n = 20
throat_radius = 35 # [m]
exit_mach = 3
gamma = 1.4

# get atmospheric properties
h = StandardAtmosphere.get_geopotential_altitude(altitude)
T, p, rho = StandardAtmosphere.get_properties(h)

# use method of characteristics to get minimum length nozzle contour
nozzle = Nozzle(throat_radius, exit_mach, gamma)
max_wall_angle, characteristics = nozzle.min_length_nozzle_moc(n, 5e-2)

# plotting resulting contour and characteristic nodes
plt.figure(1)

for node in characteristics:
    plt.plot(node.x_coordinate,node.y_coordinate, color='blue', marker='o')

for point in range(nozzle.contour.shape[1]):
    plt.plot(nozzle.contour[0,:], nozzle.contour[1,:], color='black')

plt.show()