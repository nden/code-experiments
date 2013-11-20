# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
"""
from astropy import time
from astropy import coordinates as coo
from astropy import units as u


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
        if units is not None:
            if len([units]) != num_axes:
                raise ValueError("Number of units does not match number of axes")
        self._units=units
        if axes_names is not None:
            if len([axes_]) != num_axes:
                raise ValueError("Number of axes names does not match number of axes")
        self._axes_names = axes_names
        self._reference_system = ref_system
        delf._reference_position = ref_pos

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

    reference_position :an instance of ReferencePosition

    reference_direction : (optional)
    """

    standard_time_scale = ["TT", "TDT", "ET", "TAI", "IAT",
                           "UTC", "GPS", "TDB", "TEB", "TCG",
                           "TCB", "LST"]

    def __init__(self, time_scale, reference_position, reference_direction=None, units="s", axis_name=None):
        assert time_scale in standard_time_scale, "Unrecognized time scale"
        if axis_name is None:
            axis_name = 'time'
        super("Time", numaxes=1, ref_pos=reference_position, units=[units], axes_names=axis_name)
        ##TODO? Is time_scale equivalent to reference_system?
        self._time_scale = time_scale


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
    offset_center : tuple???
    """
    standard_ref_frame = ["FK4", "FK5", "ICRS", '...']

    standard_reference_positions = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                    "TOPOCENTER", "MOON", "EMBARYCENTER", "RELOCATABLE",
                                    "UNKNOWNRefPos"]

    def __init__(self, reference_system, reference_position, obstime, equinox=None,
                 projection="", axes_names=["", ""], units=["",""], offset_center=None):


        super(SkyFrame, self).__init__(num_axes=2, ref_system=reference_system, ref_pos=reference_position,
                                       units=units, axes_names=axes_names)
        self._equinox = equinox
        self._projection = projection
        self._obstime = obstime
        self._offset_center = offset_center

    def transform_to(self, lat, lon, other):
        """
        Transform from the current reference system to other.
        """
        func = self._wrap_func(self, lat, lon, other)
        return func(lat, lon)

    def _wrap_func(self, lat, lon, other):
        def create_coordinates(lat, lon):
            coord_klass = getattr(coo, self.reference_system+"Coordinates")
            current_coordinates = coord_klass(lat, lon, unit=self.unit,
                                              obstime=self.obstime, equinox=self.equinox)
            return current_coordinates.transform_to(other)


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
    rest_wave : float or None or Quantity
        rest wavelength in units

    """
    spectral_ref_frame = ["WAVE", "FREQ", "ENER", "WAVEN",  "AWAV", "VRAD",
                          "VOPT", "ZOPT", "VELO", "BETA"]

    standard_reference_positions = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                    "TOPOCENTER", "LSR", "LSRK", "LSRD",
                                    "GALACTIC_CENTER", "MOON", "LOCAL_GROUP_CENTER",
                                    "EMBARYCENTER", "RELOCATABLE", "UNKNOWNRefPos"]

    def __init__(self, reference_system, reference_position, dateobs, units, axes_names=None,
                 rest_freq=None, rest_wave=None):
        if axes_names is None:
            axes_names = reference_system
        super(SpectralFrame, self).__init__(num_axes=1, ref_system=reference_system, ref_pos=reference_position,
                                            axes_names=axes_names, units=units)
        self._rest_freq = rest_freq
        self._rest_wave = rest_wave
        self._dateobs = dateobs

    def transfer_to(self, coordinate, other):
        u.

class CoordinateSystem(object):
    """
    A coordinate system has one or more frames.

    Parameters
    ----------
    name : str
        a user defined name
    frames : list
        list of frames [Time, Sky, Spectral]
    """
    def __init__(self, frames, name="CompositeSystem"):
        self.frames = frames

class PixelCoordinateSystem(CoordinateSystem):
    def __init__(self, reference_pixel):
        self._refpix = reference_pixel


