from math import e # Euler's constant 

class StandardAtmosphere:
    '''
    Class for calculating accurate atmospheric properties throughout 
    different regions of the Earth's atmosphere in SI units

    John Anderson Introduction to Flight (Eighth Edition)
    '''
    
    # constants
    global h0, T0, p0, rho0, a, g0, R, Re
    h0 = 0                     # sea level altitude [m]
    T0 = 288.16                # temperature at sea level [K]
    p0 = 101325                # pressure at sea level [Pa]
    rho0 = 1.2250              # density at sea level [kg/m^3]
    a = -0.0065                # first gradient region lapse rate [K/m]
    g0 = 9.8052                # gravity at sea level [m/s^2]
    R = 287                    # air gas constant [J/kgK]
    Re = 6.356766 * (10 ** 6)  # radius of Earth [m]
    
    @staticmethod
    def get_properties(h):
        """Calculates atmospheric properties by dividing atmosphere into regions

        Args:
            h (float): geopotential altitude [m]

        Returns:
            float: ambient temperature at provided altitude [K]
            float: ambient pressure at provided altitude [Pa]
            float: ambient density at provided altitude [kg/m^3]
        """

        if h == 0:                  # if sea level
            return T0, p0, rho0
        elif h < 11000:             # if under 11 km
            return StandardAtmosphere.gradient_region(h)
        elif h <= 25000:            # if between 11-25km
            T, p, rho = StandardAtmosphere.gradient_region(11000)
            dh = h - 11000          # difference in altitude from last region
            return StandardAtmosphere.isothermal_region(dh, p, T, rho)

    @staticmethod
    def get_geopotential_altitude(hg):
        """Calculates geopotential altitude from geometric altitude

        Args:
            hg (float): geometric altitude [m]

        Returns:
            float: geopotential altitude for provided geometric altitude [m]
        """
       
        if ((Re + hg) == 0):
            raise(ZeroDivisionError("get_geopotential_altitude: Division by Zero"))
        h = hg * (Re / (Re + hg))
        return h                        # where h is geopotential altitude in m

    @staticmethod
    def gradient_region(h):
        """Calculates atmospheric properties for the gradient region at an altitude

        Args:
            h (float): geopotential altitude [m]

        Returns:
            float: ambient temperature at provided altitude [K]
            float: ambient pressure at provided altitude [Pa]
            float: ambient density at provided altitude [kg/m^3]
        """
        T1 = T0
        p1 = p0
        rho1 = rho0
        h1 = h0

        T = T1 + a * (h - h1)
        p = p1 * (T / T1) ** ((-g0) / (a * R))
        rho = rho1 * (T / T1) ** (-(g0 / (a * R) + 1))

        return T, p, rho                # in SI units

    @staticmethod
    def isothermal_region(dh, p1, T1, rho1):
        """Calculates atmospheric properties for the isothermal region at an altitude 
           based on values from preceding altitude

        Args:
            dh (float): change in geopotential altitude from preceding altitude [m]
            p1 (float): pressure at preceding altitude [Pa]
            T1 (float): temperature at preceding altitude [K]
            rho1 (float): density at preceding altitude [kg/m^3]

        Returns:
            float: ambient temperature at provided altitude [K]
            float: ambient pressure at provided altitude [Pa]
            float: ambient density at provided altitude [kg/m^3]
        """

        T = T1
        p = p1 * e ** ((-1 * (g0 / (R * T))) * (dh))
        rho = rho1 * e ** ((-1 * (g0 / (R * T))) * (dh))

        return T, p, rho                # in SI unit


# validation
# for h_g in range (25): #step in km  
#     h = StandardAtmosphere.get_geopotential_altitude(h_g*1000) 
#     print(h_g*1000, round(h,0), StandardAtmosphere.get_properties(h)) 