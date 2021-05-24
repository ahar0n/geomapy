from cartography import Converter
from geodesygeometry import Ellipsoid

wgs84 = Ellipsoid(6378137, 298.257223563)

miConversor = Converter(wgs84, -69, 0.9996, 10000000, 500000)

lat, lon = miConversor.tm2geo(6958579.443, 363102.736)
print('lat: {}, lon: {}'.format(round(lat, 8), round(lon, 8)))

n, e = miConversor.geo2tm(-27.48952286, -70.38577473)
print('n: {}, e: {}'.format(round(n, 3), round(e, 3)))
