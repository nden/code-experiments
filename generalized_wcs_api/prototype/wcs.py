# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
#from astropy import wcs
#from astropy.io import fits
import pyfits as fits
from . import transforms

class ModelDimensionalityError(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message

class WCS(object):
    """
    Base WCS class

    Parameters
    ----------
    transform : astropy.modeling.Model  or a callable
        a callable which performs the transformation

    output_coordinate_system : astropy.wcs.CoordinateSystem

    """
    def __init__(self, transform, output_coordinate_system, input_coordinate_system='pixels'):
        self._transform = transform
        self._output_coordinate_system = output_coordinate_system
        self._units = self.output_coordinate_system.units
        self._input_coordinate_system = PixelCoordinateSystem()

    def select(self, label):
        if isinstance(self.transform, transforms.SelectorModel):
            return self.transform._selector[label]
        else:
            raise ModelDimensionalityError("This WCS has only one transform.")

    def __call__(self, *args):
        """
        Performs the forward transformation pix --> world.
        """
        output = self.transform(*args)
        return output

    def invert(self, args):
        try:
            inverse_transform = self.transform.inverse()
            return inverse_transform(*args)
        except NotImplementedError:
            return self.transform.invert(*args)

class FITSWCS(WCS):
    """
    Implements FITS WCS; uses `~astropy.wcs`.

    Parameters
    ----------
    fits_file : string or fits.HDUList object

    ext : int
        extension number
    """
    def __init__(self, fits_file, ext=None):
        """
        Creates an `~astropy.wcs.WCS` object and uses its `all_pix2world`
        method as a transform.
        Also constructs a `CoordSystem` object from it.

        regions mask is None in this case.
        """
        self._header, self._fobj, self._closeobj = self._get_input(fits_file, ext)
        self.fitswcsobj = wcs.WCS(self._input)
        if self._closeobj:
            self._fobj.close()
        self._transform = self.fitswcsobj.all_pix2world
        self._coordinate_system = self.create_coordinate_system()

    def __call__(self, *args):
        """
        The forward transform uses `all_pix2world` method of the object
        """
        return self._transform(args, 1)

    def create_coordinate_system(self):
        ref_system = self.fitswcsobj.wcs.radesys
        unit = self.fitswcsobj.wcs.cunit
        ctype = self.fitswcsobj.wcs.ctype
        names = [label.split('-')[0] for label in ctype]
        projcode = self.get_projcode(ctype)[0].upper()
        flavor=coordinate_systems.CoordinateFlavor(naxes=2, flavor="SPHERICAL", unit=unit,
                                   axes_names=names)
        return coordinate_systems.SpaceFrame('Sky', ref_system, coordinate_flavor=flavor,
                                              projection=projcode)

    def get_projcode(self, ctype):
        return [ctype[0][5:8], ctype[1][5:8]]

    def __invert(self, args):
        """
        The inverse transform uses the `all_sky2pix` method.
        """
        return self._fitsobj.all_world2pix(args, 1)

    def _get_input(self, filename, ext):

        closeobj = False
        if isinstance(filename, basestring):
            if is_fits(filename):
                fobj = fits.open(filename)
                closeobj = True
                header = self._get_header(fobj, ext)
            else:
                raise ValueError("Expected a FITS file.")
        elif isinstance(filename, fits.HDUList):
            fobj = filename
            header = self._get_header(fobj, ext)
        elif isinstance(filename, Header):
            header = filename
            fobj = None
        else:
            raise ValueError("Expected a FITS file name, a FITS file object or a FITS header")
        return header, fobj, closeobj

    def _get_header(self, fobj, ext):
        if ext is not None:
            assert isinstance(ext, [int, tuple]), "Expected a valid FITS extension."
        else:
            raise ValueError("'ext' parameter is required when 'fits_file' is a file"
                             " object or a filename")
        return fobj[ext]

