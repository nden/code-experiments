fomr __future__ import division
import numpy as np

class SelectorModel(Model):
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

class RegionsSelector(SelectorModel):
    """
    Parameters
    ----------
    labels : list of strings or ints
        region IDs
    transforms : list of transforms
        a transform is an instance of modeling.Model or a callable
        which performs the transformation from the region's coordinate system
        to some other coordinate system
    regions : string
        a JSON file with regions definitions
    mask_shape : tuple
        the shape of the mask

    """
    def __init__(self, labels, transforms, regions_def, mask_shape):
        super(RegionsSelector, self).__init__(labels, transforms)
        #create_regions_mask_from_json is defined in the regions API.
        self.regions_map = regions.create_regions_map_from_json(mask_shape, regions, schema)

    def __call__(self, x, y):
        """
        Transforms x and y using self._selector and self.regions_map

        Parameters
        ----------
        x, y : arrays or scalars
            x, y coordinates
        """
        result = self.regions_mask[x, y]
        unique_regions = np.unique(result)
        for i in unique_regions:
            indices = (result==i)
            transform = self._selector[i]
            result[indices]=transform(x[indices], y[indices])
        return result

    def inverse(self, x, y):
        """
        Need to invert the regions and construct an InverseSelectorModel
        """