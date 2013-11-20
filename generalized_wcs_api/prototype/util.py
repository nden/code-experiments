"""
Utility function for WCS

"""
from __future__ import division, print_function
try:
    from astropy import time
    HAS_TIME = True
except ImportError:
    HAS_TIME = False
#from astropy import coordinates

#these ctype values do not include yzLN and yzLT pairs
ctype_pairs = {"equatorial": ["RA--", "DEC-"],
               "ecliptic": ["ELON", "ELAT"],
               "galactic": ["GLON", "GLAT"],
               "helioecliptic": ["HLON", "HLAT"],
               "supergalactic": ["SLON", "SLAT"],
               "spec": specsystems
               }

radesys = ['ICRS', 'FK5', 'FK4', 'FK4-NO-E', 'GAPPT', 'GALACTIC']

def get_csystem(header):
    ctypesys = "UNKNOWN"
    ctype = [header["CTYPE1"], header["CTYPE2"]]
    radesys = header.get("RADESYS", None)
    equinox = header.get("EQUINOX", None)
    epoch = header.get("EPOCH", None)
    dateobs = header.get("MJD-OBS", header.get("DATE-OBS", None))
    cs = [ctype[0][:4], ctype[1][:4]]

    for item in ctype_pairs.items():
        if cs[0] in item[1]:
            assert cs[1] in item[1], "Inconsistent coordinate system type in CTYPE"
            ctypesys = item[0]
    if ctypesys == 'spec':
        #try to get the rest of the kw that define a spectral system from the header
        return csystems.SpectralCoordSystem(cs, **kwargs)
    if ctypesys not in ['equatorial', 'ecliptic']:
        return coordinates.__getattribute__(sky_systems_map[ctypesys])(0., 0., equinox=equinox,
                        obstime=dateobs, unit=units)
        #return ctypesys, radesys, equinox
    else:
        if radesys is None:
            if equinox is None:
                radesys = "ICRS"
            else:
                if equinox < 1984.0:
                    radesys = 'FK4'
                else:
                    radesys = 'FK5'
        if radesys in ['FK4', 'FK4-NO-E']:
            if radesys == 'FK4-NO-E':
                assert ctypesys != "ecliptic", (
                    " Inconsistent coordinate systems: 'ecliptic' and 'FK4-NO-E' ")
            if equinox is None:
                if epoch is None:
                    equinox = "1950.0"
                else:
                    equinox = epoch
        elif radesys == 'FK5':
            if equinox is None:
                if epoch is None:
                    equinox = "2000.0"
                else:
                    equinox = epoch
        elif radesys == 'GAPPT':
            assert dateobs is not None, "Either 'DATE-OBS' or 'MJD-OBS' is required"
            equinox = dateobs
    if HAS_TIME:
        if equinox is not None:
            equinox = time.Time(equinox, scale='utc')
        if dateobs is not None:
            dateobs = time.Time(dateobs, scale='utc')
    return ctypesys, radesys, equinox, dateobs
    #units = header.get('CUNIT*', ['deg', 'deg'])
    #return coordinates.__getattribute__(sky_systems_map[radesys])(0., 0., equinox=equinox,
    #                obstime=dateobs, unit=units)

def get_projcode(header):
    ctype = header.get("CTYPE*").value()
    if not ctype:
        return None
    else:
        return [ctype[0][5:8], ctype[1][5:8]]

def populate_meta(header, regions, regionsschema):
    fits_keywords = []
    if regions is None:
        fits_keywords.append('WREGIONS')
    if regionsschema is None:
        fits_keywords.append('REGSCHEM')

    meta = {}
    meta['NAXIS'] = header.get('NAXIS', 0)
    for i in range(1, meta['NAXIS']+1):
        name = 'NAXIS'+str(i)
        meta[name] = header.get(name)
        name = 'CUNIT'+str(i)
        meta[name] = header.get(name, "deg")

    for key in fits_keywords:
        meta[key] = header.get(key, None)
    return meta

specsystems = ["WAVE", "FREQ", "ENER", "WAVEN", "AWAV",
               "VRAD", "VOPT", "ZOPT", "BETA", "VELO"]

sky_systems_map = {'ICRS': 'ICRSCoordinates',
                                    'FK5': 'FK5Coordinates',
                                    'FK4': 'FK4Coordinates',
                                    'FK4NOE': 'FK4NoETermCoordinates',
                                    'GAL': 'GalacticCorrdinates',
                                    'HOR': 'HorizontalCoordinates'
                                }