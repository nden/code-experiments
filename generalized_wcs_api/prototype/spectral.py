# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and frameworks for coordinate objects.
"""

from abc import ABCMeta, abstractproperty, abstractmethod

import copy
from astropy.coordinates.transformations import *
from astropy.coordinates import transformations
from astropy.extern import six
from astropy import units as u

@six.add_metaclass(ABCMeta)
class SpectralCoordinatesBase(object):

    @abstractmethod
    def __init__(self, *args, **kwargs):
        if register:
            self.register()


    def transform_to(self, tosys):
        """
        Transform this coordinate to a new system.

        Parameters
        ----------
        tosys : class
            The system to transform this coordinate into.

        Returns
        -------
        transcoord
            A new object with this coordinate represented in the `tosys` system.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """

        #from .transformations import master_transform_graph
        #from .errors import ConvertError

        if tosys is self.__class__:
            return copy.deepcopy(self)

        trans = master_transform_graph.get_transform(self.__class__, tosys)
        #if trans is None:
        #    raise ConvertError('Cannot transform from {0} to '
        #                       '{1}'.format(self.__class__, tosys))
        return trans(self)

@transformations.coordinate_alias('wave')
class WAVE(SpectralCoordinatesBase):
    def __init__(self, sp_coord, ref_pos, unit, rest_wave=None):
        super(WAVE, self).__init__()
        if not isinstance(unit, u.Unit):
            unit = u.Unit(unit)
        if not isinstance(rest_wave, u.Quantity):
            rest_wave = rest_wave * unit
        self._rest_wave = rest_wave
        self._unit = unit
        self._ref_pos = ref_pos
        self._sp_coord = sp_coord * self._unit

def spectral_transformation(fromsys, tosys, rest_freq=None, equivalencies=None):
    return fromsys.unit.to(tosys.unit, equivalencies=equivalencies)

#def _wave_to_freq(
def optical_velocity(self, unit):
    rest_freq = self._rest_wave.to(u.Hz, equivalencies=u.spectral())
    optical_equiv = u.doppler_optical(rest_freq)
    return self._sp_coord.to(unit, equivalencies=u.optical_equiv())
