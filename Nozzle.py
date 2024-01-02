from math import tan
from MoCHelpers import MoCHelpers as moc
import numpy as np

class Nozzle:
    '''
    Class to store and calculate flow properties for a nozzle
    '''

    def __init__(self, throat_radius, exit_mach, gamma = 1.4):
        self.gamma = gamma
        self.throat_radius = throat_radius
        self.exit_mach = exit_mach
        self.contour = None

    def min_length_nozzle_moc(self, num_lines = 20, initial_flow_angle = 5e-2):
        """Calculates the minimum length contour for the nozzle object with steady, inviscid, 
        irrotational, isentropic, and supersonic flow using the Method of Characteristics (MoC)

        (all angles with respect to horizontal)

        Args:
            num_lines (int, optional): the number of characteristic lines (Defaults to 20)
            initial_flow_angle (float, optional): a small flow deflection angle to start off flow at throat
                                                  (Defaults to 5e-2)

        Returns:
            float: the maximum wall angle of the nozzle contour
            CharacteristicNode list: contains all characteristic nodes and its attributes
        """
        class CharacteristicNode:
            '''
            Class to store flow properties for an individual characteristic node

            (all angles with respect to horizontal)
            '''

            def __init__(self,index, type = "internal"):
                self.index = index               # index of node
                self.type = type                 # type of node {"internal","wall","centerline"}
                self.flow_angle = None           # flow deflection angle [deg]
                self.pm_angle = None             # Prantl-Meyer angle [deg]
                self.mach = None                 # Mach number
                self.mach_angle = None           # Mach angle [deg]
                self.right_characteristic = None # angle of right-running characteristic line (C-) 
                self.left_characteristic = None  # angle of left-running characteristic line (C+)
                self.x_coordinate = None         # x-coordinate
                self.y_coordinate = None         # y-coordinate

        # valid input checks
        if (num_lines < 3):
            raise(ValueError("min_length_nozzle_moc: Number of characteristic lines must be > 3. (current value: " +  str(num_lines) + ")"))
        elif (self.throat_radius <= 0):
            raise(ValueError("min_length_nozzle_moc: Throat radius must be > 0. (current value: " + str(self.throat_radius) + ")"))
        elif (self.exit_mach <= 1):
            raise(ValueError("min_length_nozzle_moc: Throat Mach number must be > 1. (current value: " + str(self.exit_mach) + ")"))
        elif (self.gamma <= 1):
            raise(ValueError("min_length_nozzle_moc: Specific heat ratio must be > 1. (current value: " + str(self.gamma) + ")"))

        # initialize list of characteristic nodes
        characteristics = [CharacteristicNode(0,"wall")]
        num_nodes = int(num_lines + (num_lines * ((num_lines + 1) / 2)))
        j = 1 + num_lines               # running index for when wall point is encountered
        k = 0                           # running count of which refleciton 
        
        # assign node types
        for i in range(1,num_nodes + 1):
            node = CharacteristicNode(i)
            # wall nodes
            if (i == j + k):
                node.type = "wall"
                k += 1
                j += num_lines - k
            # centerline nodes
            if (i == 1):
                node.type = "centerline"
            else:
                if (characteristics[i-1].type == "wall"): # if previous node was a wall point then next node must be a centerline point
                    node.type = "centerline"
            
            characteristics.append(node)
            
        # find max wall angle and angle intervals
        max_wall_angle = moc.prantl_meyer_function_PMangle(self.exit_mach, self.gamma,"deg") / 2
        angle_intervals = np.linspace(0,max_wall_angle,num_lines)

        # method of characteristics (MoC)

        # infinitesimally ahead of the throat (node "0")
        node = characteristics[0]
        node.flow_angle = initial_flow_angle
        node.pm_angle = node.flow_angle
        node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", 1.05)
        node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
        # node.right_characteristic = 
        # node.left_characteristic = 0
        node.x_coordinate = 0
        node.y_coordinate = self.throat_radius      
        
        # first reflection (nodes 1 through the number of characteristic lines + 1)
        for i in range(1, num_lines + 2):
            node = characteristics[i]
            adjacent_node = characteristics[0]
            previous_node = characteristics[i - 1]

            # centerline node (node 1)
            if (node.type == "centerline"):
                node.flow_angle = 0
                node.pm_angle = previous_node.flow_angle + previous_node.pm_angle
                node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", 1.05)
                node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((previous_node.flow_angle - previous_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = 0
                node.x_coordinate = self.throat_radius * tan((90 + node.right_characteristic) * (np.pi / 180))
                node.y_coordinate = 0
            # internal nodes
            elif (node.type == "internal"):
                node.flow_angle = angle_intervals[i - 1]
                node.pm_angle = node.flow_angle
                node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", previous_node.mach)
                node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
                tv_ax = (2 * angle_intervals[i - 1] + (previous_node.flow_angle - previous_node.pm_angle)) / 2
                M_ax = moc.prantl_meyer_function_Mach(tv_ax, self.gamma, "deg")
                mu_ax = moc.mach_angle_equation(M_ax, "deg")
                node.right_characteristic = ((tv_ax - mu_ax) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = moc.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
            # wall node
            elif (node.type == "wall"):
                node.flow_angle = previous_node.flow_angle
                node.pm_angle = node.flow_angle
                node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", previous_node.mach)
                node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = max_wall_angle
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = moc.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
            
            
        # remaining nodes
        j = 0   # running count of which reflection 
        k = 1   # iterator for angle_intervals
        for i in range(num_lines + 2,num_nodes + 1):
            node = characteristics[i]
            previous_node = characteristics[i - 1]
            adjacent_node = characteristics[i - (num_lines - j)]
            adjacent_node_2 = characteristics[i - (num_lines - j + 1)]

            # centerline nodes
            if (node.type == "centerline"):
                node.flow_angle = 0
                node.pm_angle = adjacent_node.flow_angle + adjacent_node.pm_angle
                node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", adjacent_node.mach)
                node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((adjacent_node.flow_angle - adjacent_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = 0
                node.x_coordinate, node.y_coordinate = moc.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,adjacent_node_2.x_coordinate,adjacent_node_2.y_coordinate, angles_unit="deg")
            # internal nodes
            elif (node.type == "internal"):
                node.flow_angle = angle_intervals[k]
                node.pm_angle = ((adjacent_node.flow_angle + adjacent_node.pm_angle) + (-previous_node.flow_angle + previous_node.pm_angle)) / 2
                node.mach = moc.prantl_meyer_function_Mach(node.pm_angle, self.gamma, "deg", previous_node.mach)
                node.mach_angle = moc.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((adjacent_node.flow_angle - adjacent_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = moc.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
                k += 1
            # wall nodes
            elif (node.type == "wall"):
                node.flow_angle = previous_node.flow_angle
                node.pm_angle = previous_node.pm_angle
                node.mach = previous_node.mach
                node.mach_angle = previous_node.mach_angle
                node.right_characteristic = previous_node.flow_angle
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = moc.xy_intersection(node.right_characteristic, adjacent_node.x_coordinate, adjacent_node.y_coordinate, node.left_characteristic, previous_node.x_coordinate, previous_node.y_coordinate, angles_unit="deg")
                k = 1
                j += 1

        # separate wall points 
        wall_x = []
        wall_y = []
        for node in characteristics:
            if (node.type == "wall"):
                wall_x.append(node.x_coordinate)
                wall_y.append(node.y_coordinate)
        
        # perform cubic spline  interpolation to fill in a smooth contour from previous wall points
        data = np.array([wall_x, wall_y], dtype=float)
        contour_x, contour_y = moc.cubic_spline_interpolation(data, num_nodes * 3)

        # assign new contour coordinates to nozzle object
        self.contour = np.array([contour_x, contour_y], dtype=float) 
        
        return max_wall_angle, characteristics


# validation
# import matplotlib.pyplot as plt

# # use method of characteristics to get minimum length nozzle contour
# nozzle = Nozzle(throat_radius=35, exit_mach=3, gamma=1.4)
# max_wall_angle, characteristics = nozzle.min_length_nozzle_moc(num_lines=20, initial_flow_angle=5e-2)

# # plotting resulting contour
# plt.figure(1)
# for node in characteristics:
#     plt.plot(node.x_coordinate,node.y_coordinate, color='blue', marker='o')
#     plt.text(node.x_coordinate, node.y_coordinate, str(node.index), fontsize=12, ha='right', va='bottom')  # Annotating the point with its number
#     print("Node " + str(node.index) + ": (x,y) = (" + str(node.x_coordinate) + ", " + str(node.y_coordinate) + ")")

# for point in range(nozzle.contour.shape[1]):
#     plt.plot(nozzle.contour[0,:], nozzle.contour[1,:])
# plt.show()