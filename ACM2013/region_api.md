Regions
-------

This file describes the region api.

Regions are necessary to describe multiple WCSs in an observation,
for example IFUs. This proposal follows the Region definitions in the
[STC document](http://www.ivoa.net/documents/PR/STC/STC-20050315.html).

* Note : This has some overlap with the aperture definitions in photutils.
Should we consider merging this work?

A Region is defined in the context of a coordinate system. Subclasses should
define `__contains__` method and a `scan` method.

    class Region(object):
        """
        Base class for regions.

        Parameters
        -------------
        rid : int or string
            region ID
        coordinate_system : astropy.wcs.CoordinateSystem instance or a string
            in the context of WCS this would be an instance of wcs.CoordinateSysem
        """
        def __init__(self, rid, coordinate_system):
            self._coordinate_system = coordinate_system
            self._rid = rid

        def __contains__(self, x, y):
            """
            Determines if a pixel is within a region.

            Parameters
            ----------
            x,y : float
                x , y values of a pixel

            Returns
            -------
                True or False

            Subclasses must define this method.
            """

        def scan(self, mask):
            """
            Sets mask values to region id for all pixels within the region.
            Subclasses must define this method.

            Parameters
            ----------
            mask : ndarray
                a byte array with the shape of the observation to be used as a mask

            Returns
            -------
            mask : array where the value of the elements is the region ID or 0 (for
                pixels which are not included in any region).
            """

An example of a Polygon region class

    class Polygon(Region):
        """
        Represents a 2D polygon region with multiple vertices

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
            self._rid = rid
            self._vertices = vertices
            self.coord_system = coord_system

        def __contains__(self, x, y):

        def scan(mask):




`create_regions_mask` is a function which creates a byte array with regions labels

    def create_regions_mask(mask_shape, regions_def, regions_schema=None):
        """
        Given a JSON file with regions definitions and a schema, this function

        - creates a byte array of shape mask_shape
        - reads in and validates the regions definitions
        - scans each region and marks each pixel in the output array with the region ID
        - returns the mask

        Given 2 regions on a detector of size 10x10 px, the mask may look like this:

        0 0 0 0 0 0 0 0 0 0
        1 1 1 1 1 1 1 1 1 1
        1 1 1 1 1 1 1 1 1 1
        1 1 1 1 1 1 1 1 1 1
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 0 0 0 0 2 2 2 2
        2 0 0 0 0 0 0 0 2 2

        """

JSON is used for regions definitions. The schema for a Set of Regions and a Polygon region is
defined below. It will be expanded with other kinds of regions.

    {"$schema": "http://json-schema.org/draft-03/schema#",
    "title": "A set of Regions",
    "type": "array",
    "items":{
    "type": "object",
    "title": "Polygon Region",
    "description": "A representation of a rectangular region on a detector",
    "properties": {
        "id": {
            "type": ["integer", "string"],
            "required": true,
            "unique": true
        },
        "coordinate_system": {
            "type": "object",
            "properties":{
                "name": {
                    "type": "string",
                    "enum": ["Cartesian", "Sky"]
                }
            },
            "required": true
        },
        "vertices": {
            "type": "array",
            "description": "Array of vertices describing the Polygon",
            "items": {
                "type": "array",
                "description": "Pairs of (x, y) coordinates, representing a vertex",
                "items": {
                    "type": "integer"
                },
                "minItems": 2,
                "maxItems": 2
            },
            "minItems": 4,
            "required": true
        }
       }
      }
    }
