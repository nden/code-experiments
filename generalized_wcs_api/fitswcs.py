from astropy.modeling import models
from astropy.modeling.core import Model
from astropy.modeling.parameters import Parameter
from astropy import units as units
from astropy import coordinates as coord
from astropy import wcs as astwcs

from gwcs import wcs
from gwcs import coordinate_frames as cf

__version__ = "2.0.0"

class FitsWcsPix2World(Model):
    """
    An attempt to integrate fits wcs with gwcs.

    Parameters
    ----------
    fitswcs : `~astropy.wcs.Wcs`
        Fits WCS object.
    origin : int
        0 or 1, the origin to use when evaluating the FITS WCS object.
        Default is 1. because of FITS.
    """
    def __init__(self, fitswcs, origin=1, **kwargs):
        self.fitswcs = fitswcs
        self.origin = origin
        self.inputs = tuple('x{0}'.format(ind) for ind in list(range(self.fitswcs.naxis)))
        self._outputs = tuple('x{0}'.format(ind) for ind in list(range(len(self.fitswcs.wcs.ctype))))
        super(FitsWcsPix2World, self).__init__(name='pix2world', **kwargs)

    def evaluate(self, x, y):
        return self.fitswcs.wcs_pix2world(x, y, self.origin)

    def inverse(self):
        return FitsWcsWorld2Pix(self.fitswcs, origin=self.origin)

class FitsWcsWorld2Pix(Model):

    def __init__(self, fitswcs, origin=1, **kwargs):
        self.fitswcs = fitswcs
        self.origin = origin
        self.inputs = tuple('x{0}'.format(ind) for ind in list(range(self.fitswcs.naxis)))
        self.outputs = tuple('x{0}'.format(ind) for ind in list(range(len(self.fitswcs.wcs.ctype))))
        super(FitsWcsWorld2Pix, self).__init__(name='world2pix', **kwargs)

    def evaluate(self, x, y):
        return self.fitswcs.wcs_world2pix(x, y, origin=self.origin)


class FITSWCS(Model):
    """

    Parameters
    ----------
    crpix : ndarray-like
        CRPIX values
    crval : ndarray-like
        CRVAL values
    ctype : ndarray-like
        CTYPE values

    """
    crpix = Parameter(default=0.)
    crval = Parameter(default=0.)
    pc = Parameter()
    ctype = Parameter()


    def __init__(crpix, crval, pc, cdelt, ctype, **kwargs):
        pass
