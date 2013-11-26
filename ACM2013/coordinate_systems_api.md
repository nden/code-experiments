Coordinate Systems
------------------

The proposed API for coordinate systems is based on the
[STC document](http://www.ivoa.net/documents/PR/STC/STC-20050315.html)
and [WCS Paper III](http://fits.gsfc.nasa.gov/fits_wcs.html)
and is a high level description of coordinate system classes.

[Prototype implementation (incomplete)] ( https://github.com/nden/code-experiments/blob/master/generalized_wcs_api/prototype/coordinate_systems.py)

Note: Some of this may be replaced (or merged) with the high level coordinate class
being worked on although it may be that it is only a space coordinate class.

A coordinate system consists of one or more coordinate frames which
describe a reference position and a reference frame.

The base class of all frames is `CoordinateFrame`. It is kept intentionally
very simple so that it can be easily extended.

    class CoordinateFrame(object):
        """
        Base class for CoordinateFrames

        Parameters
        ----------
        num_axes : int
            number of axes
        ref_system : str
            reference system (see subclasses)
        ref_pos : str or ``ReferencePosition``
            reference position or standard of rest
        units : list of units
            unit for each axis
        axes_names : list of str
            names of the axes in this frame
        name : str
            name (alias) for this frame
        """
        def __init__(self, num_axes, ref_system=None, ref_pos=None, units=None, axes_names=None, name=None):
            """ Initialize a frame"""

        def transform_to(self, other):
            """
            Transform from the current reference system to other
            """

The following basic frames are defined as subclasses of `CoordinateFrame`:

    class TimeFrame(CoordinateFrame):
        """
        Time Frame

        Parameters
        ----------
        scale : str
            time scale, one of ``standard_time_scale``
        format : str
            time representation
        obslat : float or ``astropy.coordinates.Latitude``
            observer's latitude
        obslon : float or ``astropy.coordinates.Longitude``
            observer's longitude
        units : str or ``astropy.units.Unit``
            unit for each axis
        axes_names : list of str
            names of the axes in this frame 
        """

        standard_time_scale = time.Time.SCALES

        time_format = time.Time.FORMATS

        def __init__(self, scale, format, obslat=0., obslon=0., units="s", axes_names="Time", name="Time"):

            
        def transform_to(self, val, tosys, val2=None, precision=None, in_subfmt=None,
                         out_subfmt=None, copy=False):
            """Transforms to a different scale/representation"""

    class SkyFrame(CoordinateFrame):
        """
        Sky Frame Representation
        
        Parameters
        ----------
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
        standart_ref_frame = ["FK4", "FK5", "ICRS", '...']

        def __init__(self, reference_system, reference_position, obstime, equinox=None,
                 projection="", axes_names=["", ""], units=["",""], distance=None,
                 obslat=None, obslon=None, name=None):

        def transform_to(self, lat, lon, other, distance=None):
            """
            Transform from the current reference system to other.
            """


    class SpectralFrame(CoordinateFrame):
        """
        Represents Spectral Frame

        Parameters
        ----------
        reference_frame : str
            one of spectral_ref_frame
        reference_position : ``ReferencePosition``
        obstime : time.Time instance
        units : str or units.Unit instance
            units for the spectral frame
        axes_names : str
            spectral axis name
            If None, reference_system is used
        rest_freq : float or Quantity (default None)
            rest frequency
        obslat : float or instance of `~astropy.coordinates.Latitude`
            observer's latitude
        obslon : float or instanceof `~astropy.coordinates.Longitude`
            observer's longitude
        """
        
        spectral_ref_frame = ["WAVE", "FREQ", "ENER", "WAVEN", "AWAV", "VRAD",
                          "VOPT", "ZOPT", "VELO", "BETA"]

        standard_reference_position = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                       "TOPOCENTER", "LSR", "LSRK", "LSRD",
                                       "GALACTIC_CENTER", "MOON", "LOCAL_GROUP_CENTER"]

    def __init__(self, reference_system, reference_position, units, axes_names=None,
                 obstime=None, rest_freq=None, obslat=None, obslon=None, name=None):

    def transform_to(self, coord, other):
        """ Transforms to other spectral coordinate system"""

    class CartesianFrame(CoordinateFrame):
        def __init__(self, num_axes, reference_position, units, projection=None, name=None, axes_names=None):
            ##TODO: reference position
            super(CartesianFrame, self).__init__(num_axes, ref_pos=reference_position, units=units,
                                                 axes_names=axes_names, name=name)

        def transform_to(self, other, *args):
            #use affine transform
            

    class DetectorFrame(CartesianFrame):
        def __init__(self, reference_pixel, units=[u.pixel, u.pixel], name=None, reference_position="Local",
                     axes_names=None):
            super(DetectorFrame, self).__init__(2, reference_position=reference_position, units=units, name=name, axes_names=axes_names)
            self._reference_pixel = 
            

A `CoordinateSystem` has one or more `CoordinateFrame` objects:

    class CoordinateSystem(object):
        """
        A coordinate system has one or more frames.

        Parameters
        ----------
        name : string
            a user defined name
        frames : list
            list of frames [Time, Sky, Spectral]
        """
        def __init__(self, frames, name="CompositeSystem"):
            self.frames = frames


    class ReferencePosition(object):
        """
        Encapsulates a reference position for a CoordFrame

        Table 1 in the STC document.

        """
        standart_reference_positions = ["GEOCENTER", "BARYCENTER", "HELIOCENTER",
                                        "TOPOCENTER", "LSR", "LSRK", "LSRD",
                                        "GALACTIC_CENTER", "MOON", "LOCAL_GROUP_CENTER",
                                        "EMBARYCENTER", "RELOCATABLE", "UNKNOWNRefPos"]

        def __init__(self, frame, name):
            if self._validate(frame, name):
                self._name = name
            else:
                raise ReferencePositionError

