"""
Describes a region.
Follows the STC region definitions.

This is an example of a polygon region.
"""
class Region(object):
    """
    Base class for regions

    A Region is defined in the context of a coordinate system
    """
    def __init__(self, rid, coordinate_system):
        self._coordinate_system = coordinate_system
        self._rid = rid

    def __contains__(self, x, y):
        """
        Parameters
        ----------
        x,y : float
            x , y values of a pixel

        Returns
        -------
            True or False

        Subclasses must define this method.
        """

class Polygon(Region):
    """
    Represents a 2D polygon region

    Parameters
    ----------
    rid : string
         polygon id
    vertices : list of (x,y) tuples or lists
         The list is ordered in such a way that when traversed in a
         counterclockwise direction, the enclosed area is the polygon.
         The last vertex must coincide with the first vertex, minimum
         4 vertices are needed to define a triangle
    coord_system : string
        coordinate system

    """
    def __init__(self, rid, vertices, coord_system="Cartesian"):
        pass

    def scan(self, mask):
        """
        Sets mask values to region id for all pixels which are within the region.

        Parameters
        ----------
        mask : ndarray
            a byte array with the shape of the observation to be used as a mask
        """