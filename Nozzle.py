class Nozzle:
    '''
    Class to store flow properties for a supersonic nozzle
    '''

    def __init__(self, throat_radius, exit_mach, gamma = 1.4):
        self.gamma = gamma
        self.throat_radius = throat_radius
        self.exit_mach = exit_mach


# validation
import StandardAtmosphere
import MethodOfCharacteristics
import matplotlib.pyplot as plt
import numpy as np

# parameters [SI Units]
altitude = 10000    # [m]
n = 20
throat_radius = 0.1 # [m]
exit_mach = 2

# atmospheric properties
h = StandardAtmosphere.StandardAtmosphere.get_geopotential_altitude(altitude)
T_inf, p_inf, rho_inf = StandardAtmosphere.StandardAtmosphere.get_properties(h)

# method of characteristics
nozzle = Nozzle(throat_radius, exit_mach, 1.4)
max_wall_angle, characteristics = MethodOfCharacteristics.MethodOfCharacteristics.min_length_nozzle_moc(nozzle, n, 0.19)

# plotting
plt.figure(1)
for node in characteristics:
    plt.plot(node.x_coordinate,node.y_coordinate, color='blue', marker='o')
    # plt.text(node.x_coordinate, node.y_coordinate, str(node.index), fontsize=12, ha='right', va='bottom')  # Annotating the point with its number
    # print("Node " + str(node.index) + ": (x,y) = (" + str(node.x_coordinate) + ", " + str(node.y_coordinate) + ")")
plt.show()