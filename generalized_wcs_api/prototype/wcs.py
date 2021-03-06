# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import inspect
from astropy import wcs
from astropy.io import fits
from . import transforms
from . import coordinate_systems

class ModelDimensionalityError(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message
from astropy.coordinates.transformations import TransformGraph
from .coordinate_systems import *
import networkx as nx


class WCS(nx.DiGraph):
    def __init__(self, transform, output_coordinate_system, input_coordinate_system=None, reference_pixel=None):
        super(NXWCS, self).__init__()
        if input_coordinate_system is None:
            input_coordinate_system = DetectorFrame(reference_pixel, name='Detector')
        self._input_coordinate_system = input_coordinate_system
        self._output_coordinate_system = output_coordinate_system
        self.add_edge(self._input_coordinate_system, self._output_coordinate_system, transform=transform)
        self._transform = transform

    def __call__(self, *args):
        transform = self.get_edge_data(self._input_coordinate_system, self._output_coordinate_system)['transform']
        return transform(*args)

    def select(self, label):
        if isinstance(self._transform, transforms.SelectorModel):
            return self._transform._selector[label]
        else:
            raise ModelDimensionalityError("This WCS has only one transform.")

    def add_transform(self, fromsys, tosys, transform):
        ##TODO: check that name is unique
        self.add_edge(fromsys, tosys, transform=transform)

    def invert(self, args):
        try:
            inverse_transform = self._transform.inverse()
            return inverse_transform(*args)
        except NotImplementedError:
            return self._transform.invert(*args)    

    def transform(self, fromsys, tosys, *args):
        """
        Parameters
        ----------
        fromsys : str
            name of `fromsys` coordinate system
        tosys : str
            name of `tosys` coordinate system
        """
        avilable_systems = self.coordinate_systems
        for c in available_systems:
            if c.name == fromsys:
                from_system = c
            elif c.name == tosys:
                to_system = c
            else:
                continue
        transform = self.get_edge_data(from_system, to_system)['transform']
        return transform(*args)

    @property
    def coordinate_systems(self):
        return self.nodes()


class FITSWCS(WCS):
    """
    Implements FITS WCS; uses `~astropy.wcs`.

    Parameters
    ----------
    fits_file : string or fits.HDUList object

    ext : int
        extension number
    """
    def __init__(self, fits_file, ext=None, key=""):
        """
        Creates an `~astropy.wcs.WCS` object and uses its `all_pix2world`
        method as a transform.
        Also constructs a `CoordSystem` object from it.

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
        unit = self.fitswcsobj.wcs.cunit #turn this into unit.Unit()
        ctype = self.fitswcsobj.wcs.ctype
        names = [label.split('-')[0] for label in ctype]
        projcode = self.get_projcode(ctype)[0].upper()
        obstime = self.date-obs # turn this into time.Time()
        return coordinate_systems.SpaceFrame(ref_system, 'BARICENTER', projection=projcode, axes_names=names, units=unit, name="world")

    def get_projcode(self, ctype):
        return [ctype[0][5:8], ctype[1][5:8]]

    def invert(self, args):
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

