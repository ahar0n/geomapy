from math import sqrt, acos, radians, cos, atan, pi, degrees, sin


class Ellipsoid:
    """
   El objeto Ellipsoid contiene los parámetros del elipsoide

   :param a: Semi eje mayor del elipsoide
   :param f: Factor de achatamiento de la relacion 1/f
   """

    def __init__(self, a: float, f: float):
        self.a = a
        self.f = 1/f
        self.b = self.polar_semiaxis()

    def polar_semiaxis(self):
        """

        :return: semi eje polar (b)
        """
        return self.a * (1 - self.f)

    def first_eccentricity(self):
        """

        :return: primera excentricidad (e)
        """
        return sqrt((self.a ** 2 - self.b ** 2) / self.a ** 2)

    def second_eccentricity(self):
        """

        :return: segunda excentricidad (e')
        """
        return sqrt((self.a ** 2 - self.b ** 2) / self.b ** 2)

    def polar_radius(self):
        """

        :return: radio de curvatura polar (c)
        """
        return self.a ** 2 / self.b

    def linear_eccentricity(self):
        """

        :return: excentricidad lineal (E)
        """
        return sqrt(self.a ** 2 + self.b ** 2)


class Radius:
    """
    Radios de curvaturas de la tierra (N) gran normal,
    (M) de la elipse meridiana, (P) radio medio de gauss y
    (P) del circulo paralelo

    :param ellipsoid: Objeto Ellipsoid
    :param latitude: latitud geodésicas expresada en grados sexagesimales
    """
    def __init__(self, ellipsoid, latitude):
        self.latitude = latitude
        self.ellipsoid = ellipsoid
        self.v = self.auxiliary_quantity_v()

    def auxiliary_quantity_v(self):
        """

        :return: cantidad auxiliar (V)
        :rtype: float
        """
        return sqrt(1 + self.ellipsoid.second_eccentricity() ** 2 * cos(radians(self.latitude)) ** 2)

    def curvature_in_the_meridian(self):
        """

        :return: M, radio de curvatura en el meridiano
        :rtype: float
        """
        return self.ellipsoid.polar_radius() / self.v ** 3

    def curvature_normal_section(self):
        """

        :return: N, radio de curvatura de la sección normal perpendicular al meridiano
        :rtype: float
        """
        return self.ellipsoid.polar_radius() / self.v

    def mean_radius(self):
        """

        :return: RM, radio de curvatura medio (gauss)
        :rtype: float
        """
        return self.ellipsoid.polar_radius() / self.v ** 2

    def parallel_circle(self):
        """

        :return: P, radio del circulo paralelo
        :rtype: float
        """
        return self.curvature_normal_section() * cos(radians(self.latitude))

class Converter:

    def __init__(self, ellipsoid: Ellipsoid):

        self.ellipsoid = ellipsoid

        pass

    def ecef2geo(self, x, y, z):
        """

        :return: {latitude, longitude, h}
        :rtype: dict
        """

        longitude = atan(y/x)

        e = self.ellipsoid.first_eccentricity()

        phi0 = atan(z/sqrt(x**2 + y**2) * (1+e**2/(1-e**2))) # radians

        while True:
            radius = Radius(self.ellipsoid, degrees(phi0))
            N = radius.curvature_normal_section()

            latitude = atan( z/sqrt(x**2 + y**2) * (1+ e**2*N*sin(phi0)/z) );
            if abs(degrees(latitude) - degrees(phi0)) < 0.0000000001:
                break
            
            phi0 = latitude

        
        radius = Radius(self.ellipsoid, degrees(latitude))

        if abs(degrees(latitude)) == 90:
            h = z/sin(latitude) - radius.curvature_normal_section()*(1-e**2);
        else:
            h = sqrt(x**2 + y**2)/cos(latitude) - radius.curvature_normal_section();
        
        return {'latitude': degrees(latitude), 'longitude': degrees(longitude), 'h': h}