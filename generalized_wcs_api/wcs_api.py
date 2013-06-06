"""
This describes a WCS data model.

It uses a mixin to handle multiple WCS in a file.

The alternative would be to have a special model (Mapper) which handles this.

"""

class WCS(HasMultipleTransforms):
    """
    WCS Data Model

    A WCS object performs two transformations:

    - forward : from pixel to output_coordinate_system
    - inverse : from the output_coordinate_system to pixel

    The WCS is specific for a detector and if necessary (e.g. IFU) provides
    mappings between regions on the detector and transorms.

    Parameters
    ----------
    transform : astropy.modeling.Model, a list of such objects or a callable
        a (possibly composite) model representing the pixel to world transformation
        a list of such objects in the case of multiple regions
        or a callable which performs the transformation

    output_coordinate_system : astropy.wcs.AstroCoordinateSystem

    regions_mask : a file or a byte array
        if transform is a list, this is a required parameter
        A regions mask is represented as a (compressed) byte array in the object but
        information about the regions can be stored originally in a file, e.g. a json file.
        More on regions in region.py.
    """

    def __init__(self, transform, output_coordinate_system, regions_mask=None):
        self.transform = transform
        self.output_coordinate_system = output_coordinate_system
        self.units = self.output_coordinate_system.units
        if isinstance(regions_maks, np.array):
            self.regions_mask = regions_mask
        else:
            self.regions_maks = region.create_regions_mask(regions_mask)

    def __call__(self, args):
        """
        Performs the forward transformatiion pix --> world.
        """
        if isinstance(self.transform, list) and len(self.transform) > 1:
            result = self.regions_mask[x, y]
            unique_regions = np.unique(result)
            for i in unique_regions:
                indices = result==i
                transform=self.get_forward_transform(i)
                result[indices]=transform(x[indices], y[indices])
            print('resut', result)
            return result
        else:
            return self.transform(x, y)

    def invert(self, args):
        if self.transform.has_inverse():
            inverse_transform = self.transform.inverse()
            return inverse_transform(*args)
        else:
            return self.transform.invert(*args)

class HasMultipleTransforms(object):
    """
    A mixin to map regions and transforms

    This class provides methods which map regions on the detector to transforms
    """
    def get_regions_mask(self, input):
        """
        A regions mask is represented in the WCS object as a byte array.

        Parameters
        ----------
        input : a byte array, a file or a fits extension
            the array may be compressed and stored in a fits extension
            if it's a file, a json file with a region description
        """

    def get_forward_mapping(self):
        """
        Given a list of transforms and a regions_mask
        it maps a transform with a mask.

        Assuming the regions are numbered 1,..,n
        this function maps them using something like zip(transforms, regions)

        Returns a dictionary {region_id: transform}
        """

    def get_forward_transform(self, region):
        """
        Returns the transform which corresponds to a specific region
        """
        return self.get_forward_mapping()[region]


    def to_world(self, x, y, **kwargs):
        """
        Performs the transformation from pixel to world coordinates using
        the appropriate transform depending on the region.

        x, y : arrays
            Pixel values
            In general they can be in different regions

        Returns
        -------
        xout, yout : world cordinates
        """


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
        Creates an astropy.WCS object and uses its `pix2world`
        method as a transform.
        Also constructs a `AstroCoordSystem` object from it.

        regions mask is None in this case.
        """