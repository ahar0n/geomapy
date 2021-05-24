from math import sin, cos, radians, degrees
from geodesygeometry import Radius


class MeridionalArc:
    """
    Longitud del arco de meridiano hasta e'10, en el ecuador para la latitud (phi), con (ee) segunda excentricidad y
    (c) radio polar.

    Algorithm: Arc Meridional (Conventional formulas)
    Blachut, T. J., Chrzanowski, A., & Saastamoinen, J. H. (1979).
    Urban Surveying and Mapping. New York, NY: Springer New York. pp. 18-20
    """

    def __init__(self, ellipsoid, latitude=None):
        self.ellipsoid = ellipsoid
        self.latitude = latitude

    def get_length(self):
        ee2 = self.ellipsoid.second_eccentricity() ** 2

        A0 = 1 - 3 / 4 * ee2 * (1 - 15 / 16 * ee2 * (1 - 35 / 36 * ee2 * (1 - 63 / 64 * ee2 * (1 - 99 / 100 * ee2))))

        if self.latitude is None:
            return A0

        else:
            A1 = 3 / 4 * ee2 * (1 - 25 / 16 * ee2 * (1 - 77 / 60 * ee2 * (1 - 837 / 704 * ee2 * (1 - 2123 / 1860 * ee2))))
            A2 = 5 / 8 * ee2 * (1 - 139 / 144 * ee2 * (1 - 1087 / 1112 * ee2 * (1 - 513427 / 521760 * ee2)))
            A4 = 35 / 72 * ee2 ** 2 * (1 - 125 / 64 * ee2 * (1 - 221069 / 150000 * ee2))
            A6 = 105 / 256 * ee2 ** 3 * (1 - 1179 / 400 * ee2)
            A8 = 231 / 640 * ee2 ** 4

            c = self.ellipsoid.polar_radius()
            latitude_rad = radians(self.latitude)

            # As = [A0*c, A1*c, A2, A4, A6, A8]

            return A0 * c * latitude_rad - A1 * c * sin(latitude_rad) * cos(latitude_rad) * (
                    1 + A2 * sin(latitude_rad) ** 2 +
                    A4 * sin(latitude_rad) ** 4 + A6 * sin(latitude_rad) ** 6 + A8 * sin(latitude_rad) ** 8)


class Converter:
    """
    Conversor de coordenadas gedesicas\TM

    :param ellipsoid: objeto Ellipsoid
    :param mc: longitud en el meridiano central
    :param k0: factor de escala en el meridiano central
    :param fn: falso norte
    :param fe: falso este
    """

    def __init__(self, ellipsoid, mc, k0, fn, fe):
        self.ellipsoid = ellipsoid
        self.mc = mc
        self.k0 = k0
        self.fn = fn
        self.fe = fe

    def geo2tm(self, latitude, longitude):
        """
        Coordenadas norte y este Tranverse de Mercator (TM)

        Algorithm: Geographical Coordinates into TM Coordinates
        Blachut, T. J., Chrzanowski, A., & Saastamoinen, J. H. (1979).
        Urban Surveying and Mapping. New York, NY: Springer New York. pp. 22-23

        :param latitude: latitud geodésicas en degrees
        :param longitude: longitu geodésicas en degrees
        :return: Coordenadas norte y este
        :rtype: tuple
        """
        latitude_rad = radians(latitude)
        radius = Radius(self.ellipsoid, latitude)
        P = radius.parallel_circle()
        ee2 = self.ellipsoid.second_eccentricity() ** 2

        a1 = P
        a2 = 1 / 2 * a1 * sin(latitude_rad)
        a3 = 1 / 6 * a1 * (-1 + 2 * cos(latitude_rad) ** 2 + ee2 * cos(latitude_rad) ** 4)
        a4 = 1 / 12 * a2 * (-1 + 6 * cos(latitude_rad) ** 2 + 9 * ee2 * cos(latitude_rad) ** 4 + 4 * ee2 ** 2 *
                            cos(latitude_rad) ** 6)
        a5 = 1 / 120 * a1 * (1 - 20 * cos(latitude_rad) ** 2 + (24 - 58 * ee2) * cos(latitude_rad) ** 4 + 72 * ee2 *
                             cos(latitude_rad) ** 6)
        a6 = 1 / 360 * a2 * (1 - 60 * cos(latitude_rad) ** 2 + 120 * cos(latitude_rad) ** 4)

        # cartesian coordinates
        B = MeridionalArc(self.ellipsoid, latitude).get_length()
        delta_lon = radians(longitude - self.mc)

        x = B + a2 * delta_lon ** 2 + a4 * delta_lon ** 4 + a6 * delta_lon ** 6
        y = a1 * delta_lon + a3 * delta_lon ** 3 + a5 * delta_lon ** 5

        north = self.fn + self.k0 * x
        east = self.fe + self.k0 * y

        return north, east

    def tm2geo(self, north, east):
        """
        Coordenadas curvilineas (lat, lon)

        Algorithm: Geographical Coordinates into TM Coordinates
        Blachut, T. J., Chrzanowski, A., & Saastamoinen, J. H. (1979).
        Urban Surveying and Mapping. New York, NY: Springer New York. pp. 24
        """

        x = (north - self.fn) / self.k0
        y = (east - self.fe) / self.k0

        A0 = MeridionalArc(self.ellipsoid).get_length()
        phi = x / (A0 * self.ellipsoid.polar_radius())
        B = MeridionalArc(self.ellipsoid, degrees(phi)).get_length()

        while abs(B - x) > 0.0005:  # tolerance precision in meters
            phi = phi + (x - B) / (A0 * self.ellipsoid.polar_radius())
            B = MeridionalArc(self.ellipsoid, degrees(phi)).get_length()

        radius = Radius(self.ellipsoid, degrees(phi))
        ee2 = self.ellipsoid.second_eccentricity() ** 2

        b1 = 1 / radius.parallel_circle()
        b2 = -1 / 2 * b1 ** 2 * sin(phi) * cos(phi) * (1 + ee2 * cos(phi) ** 2)
        b3 = -1 / 6 * b1 ** 3 * (2 - cos(phi) ** 2 + ee2 * cos(phi) ** 4)
        b4 = -1 / 12 * b1 ** 2 * b2 * (3 + (2 - 9 * ee2) * cos(phi) ** 2 + 10 * ee2 * cos(phi) ** 4 - 4 * ee2 ** 2 *
                                       cos(phi) ** 6)
        b5 = 1 / 120 * b1 ** 5 * (24 - 20 * cos(phi) ** 2 + (1 + 8 * ee2) * cos(phi) ** 4 - 2 * ee2 * cos(phi) ** 6)
        b6 = 1 / 360 * b1 ** 4 * b2 * (45 + 16 * cos(phi) ** 4)

        latitude = degrees(phi + b2 * y ** 2 + b4 * y ** 4 + b6 * y ** 6)
        longitude = degrees(radians(self.mc) + b1 * y + b3 * y ** 3 + b5 * y ** 5)

        return latitude, longitude
