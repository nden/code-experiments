# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
"""
from astropy import time
from astropy import coordinates as coo
from astropy import units as u
from astropy import utils as astu


__all__ = ['PixelCoordinateSystem', 'SkyFrame', 'SpectralFrame', 'TimeFrame',
           'CoordinateSystem']


class CoordinateFrame(object):
    """
    Base class for CoordinateFrames

    Parameters
    ----------
    system : str
        type of the frame
    num_axes : int
    axes_names : list of str
    units : list of units
    """
    def __init__(self, num_axes, ref_system=None, ref_pos=None, units=None, axes_names=None):
        """ Initialize a frame"""
        if units is not None and astu.isiterable(units):
            if len(units) != num_axes:
                raise ValueError("Number of units does not match number of axes")
        self._units=units
        if axes_names is not None and astu.isiterable(axes_names):
            if len(axes_names) != num_axes:
                raise ValueError("Number of axes names does not match number of axes")
        self._axes_names = axes_names
        self._reference_system = ref_system
        self._reference_position = ref_pos

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, val):
        self._units = val

    @property
    def axes_names(self):
        return self._axes_names

    @units.setter
    def axes_names(self, val):
        self._axes_names = val

    def transform_to(self, other):
        """
        Transform from the current reference system to other if
        the system attribute of the two matches
        """
        raise NotImplementedError("Subclasses should implement this")


    def _validate_reference_position(self):
        assert self._reference_position in self.standard_reference_position, "Unrecognized reference position"

class TimeFrame(CoordinateFrame):
    """
    Time Frame

    Parameters
    ----------
    time_scale : time scale

    scale : str
        time scale, one of ``standard_time_scale``
    format : str
        time representation

    reference_direction : (optional)
    """

    standard_time_scale = time.Time.SCALES

    time_format = time.Time.FORMATS

    def __init__(self, scale, format, obslat=0., obslon=0., units="s", axes_names="Time"):
        assert scale in self.standard_time_scale, "Unrecognized time scale"
        assert format in self.time_format, "Unrecognized time representation"
        super(TimeFrame, self).__init__(1, units=[units], axes_names=[axes_names])
        self._scale = scale
        self._format = format
        self._obslat = obslat
        self._obslon = obslon

    def transform_to(self, val, tosys, val2=None, precision=None, in_subfmt=None,
                     out_subfmt=None, copy=False):
        """Transforms to a different scale/representation"""

        t = time.Time(val, val2, format=self._format, scale=self._scale,
                      precision=precision, in_subfmt=in_subfmt, out_subfmt=out_subfmt,
                      lat=self._obslat, lon=self._obslon)
        return getattr(t, tosys)


class SkyFrame(CoordinateFrame):
    """
    Space Frame Representation

    reference_system : str
        one of standard_ref_frame
    reference_position : str, ReferencePosition instance
    obstime : time.Time instance
        observation time
    equinox : time.Time
        equinox
    ##TODO? List of supported projections per frame
    projection : str
        projection type
    axes_names : iterable or str
        axes labels
    units : str or units.Unit instance or iterable of those
        units on axes
    distance : Quantity
        distance from origin
        see ``~astropy.coordinates.distance`
    obslat : float or ~astropy.coordinates.Latitute`
        observer's latitude
    obslon : float or ~astropy.coordinates.Longitude
        observer's longitude

    """
    standard_ref_frame = ["FK4", "FK5", "ICRS", '...']

    standard_reference_positions = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                    "TOPOCENTER", "MOON", "EMBARYCENTER", "RELOCATABLE",
                                    "UNKNOWNRefPos"]

    def __init__(self, reference_system, reference_position, obstime, equinox=None,
                 projection="", axes_names=["", ""], units=["",""], distance=None,
                 obslat=None, obslon=None):


        super(SkyFrame, self).__init__(num_axes=2, ref_system=reference_system, ref_pos=reference_position,
                                       units=units, axes_names=axes_names)
        self._equinox = equinox
        self._projection = projection
        self._obstime = obstime
        self._distance = distance
        self._oblat = obslat
        self._obslon = obslon

    @property
    def obstime(self):
        return self._obstime

    @property
    def equinox(self):
        return self._equinox

    def transform_to(self, lat, lon, other):
        """
        Transform from the current reference system to other.
        """
        func = self._wrap_func(lat, lon, other)
        return func(lat, lon)

    def _wrap_func(self, lat, lon, other):
        def create_coordinates(lat, lon):
            coord_klass = getattr(coo, self._reference_system)
            current_coordinates = coord_klass(lat, lon, unit=self.units,
                                              obstime=self.obstime, equinox=self.equinox)
            return current_coordinates.transform_to(other)
        return create_coordinates


class SpectralFrame(CoordinateFrame):
    """
    Represents Spectral Frame

    Parameters
    ----------
    reference_frame : str
        one of spectral_ref_frame
    reference_position : an instance of ReferencePosition
    dateobs : time.Time instance
    units : str or units.Unit instance
        units for the spectral frame
    axes_names : str
    rest_freq : float?? or None or Quantity
        rest frequency in units
    obslat : float or instance of `~astropy.coordinates.Latitude`
        observer's latitude
    obslon : float or instanceof `~astropy.coordinates.Longitude`
        observer's longitude
    """
    spectral_ref_frame = ["WAVE", "FREQ", "ENER", "WAVEN",  "AWAV", "VRAD",
                          "VOPT", "ZOPT", "VELO", "BETA"]

    standard_reference_position = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                   "TOPOCENTER", "LSR", "LSRK", "LSRD",
                                   "GALACTIC_CENTER", "MOON", "LOCAL_GROUP_CENTER"]

    ##TODO: observer location for topocentric
    def __init__(self, reference_system, reference_position, units, axes_names=None,
                 obstime=None, rest_freq=None, obslat=None, obslon=None):

        if axes_names is None:
            axes_names = [reference_system]
        super(SpectralFrame, self).__init__(num_axes=1, ref_system=reference_system, ref_pos=reference_position,
                                            axes_names=axes_names, units=units)
        self._validate_reference_position()
        self._rest_freq = rest_freq
        self._obstime = obstime

    def rest_wavelength(self):
        return self._rest_frequency.to(u.m, equivalencies=u.spectral())

    def transfer_to(self, coordinate, other):
        return u.Quantity(coordinate).to(other.units, equivalencies=u.spectral())

class CoordinateSystem(object):
    """
    A coordinate system has one or more frames.

    Parameters
    ----------
    name : str
        a user defined name
    frames : list
        list of frames
        one of TimeFrame, SkyFrame, SpectralFrame, CoordinateFrame, CoordinateSystem
    """
    def __init__(self, frames, name="CompositeSystem"):
        self.frames = frames

class PixelCoordinateSystem(CoordinateSystem):
    def __init__(self, reference_pixel):
        frame = CoordinateFrame(2, ref_pos=reference_pixel, axes_names=['X', 'Y'])
        super(PixelCoordinateSystem, self).__init__(frames=frame, name='Pixel')

class ReferencePosition(object):

    standard_reference_position = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                   "TOPOCENTER", "LSR", "LSRK", "LSRD",
                                   "GALACTIC_CENTER", "MOON", "EMBARYCENTER",
                                   "RELOCATABLE", "UNKNOWNRefPos"]

    def __init__(self, ref_pos, lat=None, lon=None, alt=None):
        if not ref_pos in standard_reference_position:
            raise ValueError("Unrecognized reference position")

        if ref_pos == 'TOPOCENTRIC' and (lat is None or lon is None or alt is None):
            raise ValueError("Need lat, lon, alt to define Topocentric reference position")

        self._ref_pos = ref_pos

    @property
    def ref_pos(self):
        return self._ref_pos
