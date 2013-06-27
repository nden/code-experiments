class CoordinateFrame(object):
    """
    Base class for CoordinateFrames

    Parameters
    ----------
    system : string
        type of the frame
    num_axes : int
    axes_names : list of strings
    units : list of units
    """
    def __init__(self, system, num_axes, axes_names=None, units=None):
        """ Initialize a frame"""

    def transform_to(self, other):
        """
        Transform from the current reference system to other if
        the system attribute of the two matches
        """

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

    def __init__(self, time_scale, reference_position, reference_direction=None, units="s"):
        assert time_scale in standard_time_scale, "Unrecognized time scale"

        super("Time", numaxes=1, axes_names=['time'], units=[units])
        self._time_scale =time_scale
        self._reference_position = reference_position

class SkyFrame(CoordinateFrame):
    """
    Space Frame Representation
    """
    standart_ref_frame = ["FK4", "FK5", "ICRS", '...']

    def __init__(self, reference_frame, reference_position, obstime, equinox=None,
                 projection="", axes_names=["", ""], units=["",""], offset_center=None):

        super(SkyFrame, self).__init__('Space', num_axes=2, axes_names=axes_names, units=units)
        self._equinox = equinox
        self._projection = projection
        self._offset_center = offset_center
        self._reference_frame = reference_frame
        self._reference_position = reference_position

class SpectralFrame(CoordinateFrame):
    """
    Represents Spectral Frame

    Parameters
    ----------
    reference_frame : string
        one of spectral_ref_frame
    reference_position : an instance of ReferencePosition

    """
    spectral_ref_frame = ["WAVE", "FREQ", "ENER", "WAVEN",  "AWAV", "VRAD",
                        "VOPT", "ZOPT", "VELO", "BETA"]
    def __init__(self, reference_frame, reference_pos, date_obs, rest_freq=None, rest_wave=None,
                          units=""):

        super(SpectralFrame, self).__init__('Spec', num_axes=1, axes_names=[reference_frame], units=[""])
        self._reference_frame = reference_frame
        self._reference_position = reference_position
        self.rest_freq = rest_freq
        self.rest_wave = rest_wave
        self.dateobs = dateobs

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

    def _validate_ref_pos(self, frame):
        """
        validates that the reference position is allowed for the frame
        """
        return True
