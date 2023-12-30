from math import atan, sqrt, asin, tan
import numpy as np
from scipy.optimize import fsolve


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

class MethodOfCharacteristics:
    '''
    Class for calculating flow characteristics throughout 
    a supersonic process using the Method of Characteristics (MOC)

    John Anderson Fundamentals of Aerodynamics (Sixth Edition)
    ''' 

    @staticmethod
    def min_length_nozzle_moc(nozzle, num_lines = 20, initial_flow_angle = 5e-2):

        # validity checks
        if (num_lines < 3):
            raise(ValueError("min_length_nozzle_moc: Number of characteristic lines must be > 3. (current value: " +  str(num_lines) + ")"))
        elif (nozzle.throat_radius <= 0):
            raise(ValueError("min_length_nozzle_moc: Throat radius must be > 0. (current value: " + str(nozzle.throat_radius) + ")"))
        elif (nozzle.exit_mach <= 1):
            raise(ValueError("min_length_nozzle_moc: Throat Mach number must be > 1. (current value: " + str(nozzle.exit_mach) + ")"))
        elif (nozzle.gamma <= 1):
            raise(ValueError("min_length_nozzle_moc: Specific heat ratio must be > 1. (current value: " + str(nozzle.gamma) + ")"))

        # initialize list of characteristic nodes
        characteristics = [CharacteristicNode(0,"wall")]
        num_nodes = int(num_lines + (num_lines * ((num_lines + 1) / 2)))
        j = 1 + num_lines               # running index for when wall point is encountered
        k = 0                           # running count of which refleciton 
        
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
            
            # push back node
            characteristics.append(node)
            
        # find max wall angle and angle intervals
        max_wall_angle = MethodOfCharacteristics.prantl_meyer_function_PMangle(nozzle.exit_mach, nozzle.gamma,"deg") / 2
        angle_intervals = np.linspace(0,max_wall_angle,num_lines)

        # method of characteristics (moc)

        # infinitesimally ahead of the throat (node "0")
        node = characteristics[0]
        node.flow_angle = initial_flow_angle
        node.pm_angle = node.flow_angle
        node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", 1.05)
        node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
        # node.right_characteristic = 
        # node.left_characteristic = 0
        node.x_coordinate = 0
        node.y_coordinate = nozzle.throat_radius      
        
        # first reflection (nodes 1 through the number of characteristic lines + 1)
        for i in range(1, num_lines + 2):
            node = characteristics[i]
            adjacent_node = characteristics[0]
            previous_node = characteristics[i - 1]

            # centerline node (node 1)
            if (node.type == "centerline"):
                node.flow_angle = 0
                node.pm_angle = previous_node.flow_angle + previous_node.pm_angle
                node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", 1.05)
                node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((previous_node.flow_angle - previous_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = 0
                node.x_coordinate = nozzle.throat_radius * tan((90 + node.right_characteristic) * (np.pi / 180))
                node.y_coordinate = 0
            # internal nodes
            elif (node.type == "internal"):
                node.flow_angle = angle_intervals[i - 1]
                node.pm_angle = node.flow_angle
                node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", previous_node.mach)
                node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
                tv_ax = (2 * angle_intervals[i - 1] + (previous_node.flow_angle - previous_node.pm_angle)) / 2
                M_ax = MethodOfCharacteristics.prantl_meyer_function_Mach(tv_ax, nozzle.gamma, "deg")
                mu_ax = MethodOfCharacteristics.mach_angle_equation(M_ax, "deg")
                node.right_characteristic = ((tv_ax - mu_ax) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = MethodOfCharacteristics.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
            # wall node
            elif (node.type == "wall"):
                node.flow_angle = previous_node.flow_angle
                node.pm_angle = node.flow_angle
                node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", previous_node.mach)
                node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = max_wall_angle
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = MethodOfCharacteristics.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
            
            
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
                node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", adjacent_node.mach)
                node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((adjacent_node.flow_angle - adjacent_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = 0
                node.x_coordinate, node.y_coordinate = MethodOfCharacteristics.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,adjacent_node_2.x_coordinate,adjacent_node_2.y_coordinate, angles_unit="deg")
            # internal nodes
            elif (node.type == "internal"):
                node.flow_angle = angle_intervals[k]
                node.pm_angle = ((adjacent_node.flow_angle + adjacent_node.pm_angle) + (-previous_node.flow_angle + previous_node.pm_angle)) / 2
                node.mach = MethodOfCharacteristics.prantl_meyer_function_Mach(node.pm_angle, nozzle.gamma, "deg", previous_node.mach)
                node.mach_angle = MethodOfCharacteristics.mach_angle_equation(node.mach, "deg")
                node.right_characteristic = ((adjacent_node.flow_angle - adjacent_node.mach_angle) + (node.flow_angle - node.mach_angle)) / 2
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = MethodOfCharacteristics.xy_intersection(node.right_characteristic,adjacent_node.x_coordinate,adjacent_node.y_coordinate,node.left_characteristic,previous_node.x_coordinate,previous_node.y_coordinate, angles_unit="deg")
                k += 1
            # wall nodes
            elif (node.type == "wall"):
                node.flow_angle = previous_node.flow_angle
                node.pm_angle = previous_node.pm_angle
                node.mach = previous_node.mach
                node.mach_angle = previous_node.mach_angle
                node.right_characteristic = previous_node.flow_angle
                node.left_characteristic = ((previous_node.flow_angle + previous_node.mach_angle) + (node.flow_angle + node.mach_angle)) / 2
                node.x_coordinate, node.y_coordinate = MethodOfCharacteristics.xy_intersection(node.right_characteristic, adjacent_node.x_coordinate, adjacent_node.y_coordinate, node.left_characteristic, previous_node.x_coordinate, previous_node.y_coordinate, angles_unit="deg")
                k = 1
                j += 1
        
        return max_wall_angle, characteristics

    @staticmethod
    def prantl_meyer_function_PMangle(M, gamma, angles_unit = "rad"):
        """Calculates the Prantl-Meyer angle at a Mach number

        Args:
            M (float): Mach number
            gamma (float): specific heat ratio
            angles_unit (string, optional): specifies whether provided angle is in radians or degrees 
                                            {"rad","deg"} (Defaults to "rad")

        Returns:
            float: Prantl-Meyer angle [rad or deg]
        """

        A = sqrt((gamma + 1) / (gamma - 1))
        B = sqrt(M ** 2 - 1)
        PM_angle = A * atan((1 / A) * B) - atan(B)

        if (angles_unit == "rad"):
            return PM_angle
        elif (angles_unit == "deg"):
            return PM_angle * (180 / np.pi)
        else:
            raise(ValueError("prantl_meyer_function_PMangle: Unit type is not valid."))
    
    @staticmethod
    def prantl_meyer_function_Mach(PM_angle, gamma, angles_unit = "rad", M_guess = 1.2):
        """Calculates the Mach number at a Prantl-Meyer angle

        Args:
            PM_angle (float): Prantl-Meyer angle
            gamma (float): specific heat ratio
            angles_unit (string, optional): specifies whether provided angle is in radians or degrees 
                                            {"rad","deg"} (Defaults to "rad")
            M_guess (float, optional): Mach number guess (Defaults to 1.2)

        Returns:
            float: Prantl-Meyer angle [rad or deg] 
        """

        if (angles_unit == "deg"):
            PM_angle *= (np.pi / 180)
        elif (angles_unit != "rad"):
            raise(ValueError("prantl_meyer_function_Mach: Unit type is not valid."))

        A = sqrt((gamma + 1) / (gamma - 1))

        return fsolve(lambda M: PM_angle - (sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - atan(sqrt(M ** 2 - 1))), M_guess)[0]

    @staticmethod
    def mach_angle_equation(M, angles_unit = "rad"):
        """Calculates the Mach angle at a Mach number

        Args:
            M (float): Mach number
            angles_unit (string, optional): specifies whether provided angle is in radians or degrees 
                                            {"rad","deg"} (Defaults to "rad")

        Returns:
            float: Mach angle [rad or deg] 
        """

        mach_angle = asin(1 / M)

        if (angles_unit == "rad"):
            return mach_angle
        elif (angles_unit == "deg"):
            return mach_angle * (180 / np.pi)
        else:
            raise(ValueError("mach_angle_equation: Unit type is not valid."))
        
    @staticmethod
    def xy_intersection(right_characteristic, right_x, right_y, left_characteristic, left_x, left_y, angles_unit = "rad"):
        """Calculates the xy-coordinates of a point downstream where the right-running
        and left-running characteristic lines intersect using rearranged forms of the 
        point-slope formula

        Args:
            right_characteristic (float): slope of right-running characteristic line
            right_x (float): x-coordinate of point where right-running characteristic line originates from
            right_y (float): y-coordinate of point where right-running characteristic line originates from
            left_characteristic (float): slope of left-running characteristic line
            left_x (float): x-coordinate of point where right-running characteristic line originates from
            left_y (float): x-coordinate of point where right-running characteristic line originates from
            angles_unit (string, optional): specifies whether provided angle is in radians or degrees 
                                            {"rad","deg"} (Defaults to "rad")
        Returns:
            x: x-coordinate of point where right-running and left-running characteristic lines intersect
            y: y-coordinate of point where right-running and left-running characteristic lines intersect
        """

        if (angles_unit == "deg"):
            right_characteristic *= (np.pi / 180)
            left_characteristic *= (np.pi / 180)
        elif (angles_unit != "rad"):
            raise(ValueError("xy_intersection: Unit type is not valid."))

        if (tan(right_characteristic) - tan(left_characteristic) == 0):
            raise(ZeroDivisionError("xy_intersection: Division by Zero"))
        x = (right_x * tan(right_characteristic) - left_x * tan(left_characteristic) + left_y - right_y) / (tan(right_characteristic) - tan(left_characteristic))
        y = (x - right_x) * tan(right_characteristic) + right_y

        return x, y


