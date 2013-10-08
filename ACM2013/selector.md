SelectorModel
-------------

The `SelectorModel` provides mappings between a transform and some other quantity.
The obvious case is an IFU detector where each region maps to a different transform.
However, the goal is for this to be general and allow for other use cases, for example
multiple orders, slitless spectroscopy.

For this reason the base class defines only the mapping and additional functionality
is left for specific classes.

The `__call__` method is defined by subclasses and implements the real work of
figuring out how to match input coordinates with a transform.

    class SelectorModel(object):
        """
        Parameters
        ----------
        labels : a list of strings or objects
        transforms : a list of transforms
            transforms match labels
        """
        def __init__(self, labels, transforms):
            self._selector = dict(labels, transforms)

        def __call__(self, label):
           raise NotImplementedError

An example of a selector of regions of an IFU is below. It has an additional attribute
`regions_map` which is a byte array with region labels. The actual definitions of the
regions are read in from a JSON file, see the [Region API](https://github.com/nden/astropy-api/blob/generalized_wcs/generalized_wcs/region_api.md)

    class RegionsSelector(SelectorModel):
        def __init__(self, labels, transforms, regions):
            super(RegionsSelector, self).__init__(labels, transforms)
            #create_regions_mask_from_json is defined in the regions API.
            self.regions_map = regions.create_regions_map_from_json(mask_shape, regions, schema)

        def __call__(self, x, y):
            """
            Transforms x and y using self._selector and self.regions_map

            """