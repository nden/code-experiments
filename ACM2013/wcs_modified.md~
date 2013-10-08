WCS data model
--------------

Written by Nadia Dencheva (@nden) with ideas and feedback from
Perry Greenfield (@perrygrienfield) and Mike Droettboom (@mdboom) and taking into account
discussions at the latest astropy coordination meeting (Baltimore, 2012).

This document describes the WCS data model. WCS serialization is not discussed here.

Motivation for this work was discussed [previously](https://mail.google.com/mail/u/1/?ui=2&shva=1#search/perrygreenfield/13efbb0f6b4b3932).

The basic requirements for the data model are:

* Perform the transformation from pixels to world coordinates and back.
* Handle multiple transforms ( IFU detectors but not limited to that).

High level  interface
---------------------

Create a WCS object:

    wcsobj = WCS(coordinate_system, transform)

Transform pixel coordinates to world coordinates:

    ra, dec = wcsobj(x, y)

Transform from world coordinates to pixels:

    x, y = wcsobj.invert(ra, dec)

`transform` may be an instance of `SelectorModel` which maps transforms to other quantities,
for example regions on detector but not limited to this. In this case:

    wcsobj(x,y)

should be able to transform the cordinates using the correct transform.
This is handled by the `SelectorModel` class.

    wcs_for_region_3 = wcsobj.select(label=3)

returns the WCS object for region 3.

If `transform` is an instance of `SelectorModel`, the `invert` method should map to the `invert()`
methods of all individual transforms.

It should be possible to transform the output coordinates to a different coordinate system. The `~astropy.coordinates` package
is be used where possible.

    wcsobj.coordinate_system.to(ra, dec, other_system)

The special case of FITS WCS is handled through the `FITSWCS` class which inherits from `WCS`
and uses `pywcs`  to create the coordinate system object and the transforms.

    fwcs = FITSWCS(fits_file, ext, key=' ')

In more detail
--------------

* [WCS API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/wcs_api.md#wcs-api)
* [Selector Object](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/selector.md)
* [Regions API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/region_api.md)
* [Coordinate Systems API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/coordinate_systems_api.md)

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
            return self.transform(x, y)

        def invert(self, args):
            try:
                inverse_transform = self.transform.inverse()
                return inverse_transform(*args)
            except NotImplementedError:
                return self.transform.invert(*args)

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

