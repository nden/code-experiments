WCS data model
--------------

High level  interface
---------------------

Create a WCS object:

    wcsobj = WCS(transform, output_coordinate_system)

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

It is possible to transform the output coordinates to a different coordinate system using other packages in astropy

    wcsobj.coordinate_system.to(ra, dec, other_system)

The special case of FITS WCS is handled through the `FITSWCS` class which inherits from `WCS`
and uses the current `astropy.wcs`  to create the coordinate system object and the transforms.

    fwcs = FITSWCS(fits_file, ext, key=' ')

In more detail
--------------

* [WCS API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/wcs_api.md#wcs-api)
* [Selector Object](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/selector.md)
* [Regions API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/region_api.md)
* [Coordinate Systems API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/coordinate_systems_api.md)

WCS API
-------

The WCS class is represented as a directed graph whose nodes are coordinate systems and edges are transfroms
between them. This allows one to specify multiple coordinate systems and the transformations between them. 
(The [prototype](https://github.com/nden/code-experiments/blob/master/generalized_wcs_api/prototype/wcs.py) uses `networkx`.)

    class WCS(nx.DiGraph):
        """
        WCS Data Model

        Parameters
        ----------
        transform : astropy.modeling.Model  or a callable
            a callable which performs the transformation

        output_coordinate_system : astropy.wcs.CoordinateSystem
            output system  
        input_coordinate_system : astropy.wcs.CoordinateSystem
            input coordinate system
            default is `wcs.DetectorFrame`
        """

        def __init__(self, transform, output_coordinate_system, input_coordinate_system=None):
            self.transform = transform
            self.output_coordinate_system = output_coordinate_system
            self.units = self.output_coordinate_system.units

        def select(self, label):
            """
            If the transform is an instance of `SelectorModel` pick one of the available transforms
            
            Parameters
            ----------
            label : int or str
                
            """
            if isinstance(self.transform, SelectorModel):
                return self.transform[label]
            else:
                raise Error

        def add_transform(self, fromsys, tosys, transform):
            """
            Add new coordinate systems and transforms between them to this WCS
            """
            ##TODO: check that name is unique
            self.add_edge(fromsys, tosys, transform=transform)


        def transform(self, fromsys, tosys, *args):
            """
            Perform a transformation between any two systems registered with this WCS.
            
            Parameters
            ----------
            fromsys : str
                name of `fromsys` coordinate system
            tosys : str
                name of `tosys` coordinate system
            """
        
        @property
        def coordinate_systems(self):
            return self.nodes()
            
            
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

        Uses `astropy.wcs`

        Parameters
        ----------
        fits_file : string or fits.HDUList object
        ext : int
            extension number
        key : str
            alternate WCS key
            
        """
        def __init__(self, fits_file, ext=0, key=""):
            """
            Creates an `~astropy.wcs.WCS` object and uses its `all_pix2world`
            method as a transform.
            Also constructs a `CoordSystem` object from it.
            If available, adds additional coordinate systems and transformations between them 
            based on `wcs.pix2foc` and similar methods.

            """

        def __call__(self, args):
            """
            The forward transform uses `all_pix2world` method of the object
            """

        def invert(self, args):
            """
            The inverse transform uses the `all_sky2pix` method.
            """

