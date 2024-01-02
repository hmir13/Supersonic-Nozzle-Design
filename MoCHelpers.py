from math import atan, sqrt, asin, tan
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline

class MoCHelpers:
    '''
    Class containing helper functions for Method of Characteristics for aerodynamic flows

    John Anderson Fundamentals of Aerodynamics (Sixth Edition)
    ''' 

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

    @staticmethod
    def lagrange_interpolation(data, num_intervals = 1000):
        """Performs Lagrange interpolation based on given data points

        Args:
            data (numpy.ndarray): A 2D numpy array where data[0, :] contains x-coordinates 
                                and data[1, :] contains corresponding y-coordinates of the 
                                data points for interpolation
            num_intervals (int, optional): Number of intervals to use for generating the interpolated values
                                            (Defaults to 1000)

        Returns:
            numpy.ndarray: An array of x-values generated for interpolation
            numpy.ndarray: An array of y-values representing the interpolated function 
            corresponding to the x-values
        """
        n = data.shape[1]       # polynomial degree
        x = np.linspace(min(data[0]), max(data[0]), num_intervals) # domain for interpolated data
        L = np.zeros_like(x)    # interpolated values

        for i in range(n):
            product = np.ones_like(x)

            for j in range(n):
                if j != i:
                    product *= (x - data[0,j]) / (data[0,i] - data[0,j])

            L += (product * data[1,i])
        
        return x, L

    @staticmethod
    def cubic_spline_interpolation(data, num_intervals=1000):
        """
        Perform cubic spline interpolation based on given data points

        Args:
            data (numpy.ndarray): A 2D numpy array where data[0, :] contains x-coordinates 
                                and data[1, :] contains corresponding y-coordinates of the 
                                data points for interpolation
            num_intervals (int, optional): Number of intervals to use for generating the interpolated values
                                        (Defaults to 1000)

        Returns:
            numpy.ndarray: An array of x-values generated for interpolation
            numpy.ndarray: An array of y-values representing the interpolated function 
                        corresponding to the x-values
        """
        
        x_data, y_data = data[0], data[1]
        
        # Create the cubic spline interpolation function
        cs = CubicSpline(x_data, y_data)
        
        # Generate x-values based on the minimum and maximum x-values from the data
        x = np.linspace(min(x_data), max(x_data), num_intervals)
        
        # Evaluate the cubic spline at the generated x-values
        interpolated_values = cs(x)
        
        return x, interpolated_values

