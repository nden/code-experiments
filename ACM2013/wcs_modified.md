WCS data model
--------------


WCS API
-------

    class WCS(object):
        """
        WCS Data Model

        Parameters
        ----------
        transform : astropy.modeling.Model  or a callable
            a callable which performs the transformation

        output_coordinate_system : astropy.wcs.CoordinateSystem
        """

        def __init__(self, transform, output_coordinate_system):
            self.transform = transform
            self.output_coordinate_system = output_coordinate_system
            self.units = self.output_coordinate_system.units

        def select(self, label):
            if isinstance(self.transform, SelectorModel):
                return self.transform[label]
            else:
                raise Error

        def __call__(self, args):
            """
            Performs the forward transformation pix --> world.
            """
            return self.transform(*args)

        def invert(self, args):
            try:
                inverse_transform = self.transform.inverse()
                return inverse_transform(*args)
            except NotImplementedError:
                return self.transform.invert(*args)

        def add_transform(from_coord_system, to_coordinate_system, transform):
            """
            Adds a coordinate system and a transform to the graph.
            from_coord_system must already exist
            """

        def transform(from_coord_sys, to_coord_sys, *args):
            """
            Check if transform between two coordinate systems exists
            and return the transformed coordinates.
            """

FITS WCS uses `~astropy.wcs`. This takes advantage of knowledge about FITS
specific keywords and parsing the FITS header.

    class FITSWCS(WCS):
        """
        Implements FITS WCS

        Uses pywcs

        Parameters
        ----------
        fits_file : string or fits.HDUList object

        ext : int
            extension number
        """
        def __init__(self, fits_file, ext=0):
            """
            Creates an `~astropy.wcs.WCS` object and uses its `all_pix2world`
            method as a transform.
            Also constructs a `CoordSystem` object from it.

            regions mask is None in this case.
            """

        def __call__(self, args):
            """
            The forward transform uses `all_pix2world` method of the object
            """

        def __invert(self, args):
            """
            The inverse transform uses the `all_sky2pix` method.
            """

