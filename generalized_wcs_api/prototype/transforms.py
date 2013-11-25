# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import copy
import numpy as np
from astropy.modeling import *
from astropy import utils as astu
from . import region
import json

class ImagingWCS(object):
    def __init__(self, wcs_info, dist_info):
        wcsaxes = wcs_info['WCSAXES']
        crpix = wcs_info['CRPIX']
        pc = np.array(wcs_info['PC'], dtype=np.float)
        pc.shape = (2, 2)
        crval = wcs_info['CRVAL']
        phip = crval[0]
        lonp = crval[1]
        thetap = 180 # write "def compute_lonpole(projcode, l)"
        ctype = wcs_info['CTYPE']
        projcode = get_projcode(ctype)[0].upper()
        if projcode not in projections.projcodes:
            raise ValueError('Projection code %s, not recognized' %projcode)
        # create a projection transform
        self.projection = create_projection_transform(projcode)
        # Create a RotateNative2Celestial transform using the Euler angles
        self.n2c = models.RotateNative2Celestial(phip, thetap, lonp)

        # create shift transforms in x and y
        self.offx = models.ShiftModel(-crpix[0])
        self.offy = models.ShiftModel(-crpix[1])

        # create a CD/PC transform
        self.pcrot = models.MatrixRotation2D(rotmat=pc)
        if dist_info is not None:
            self.distortion = self.create_distortion_transform(dist_info)
        else:
            self.distortion = None

    def create_distortion_transform(self, distortion):
        dist_info = copy.deepcopy(distortion)
        try:
            dist_info.pop('title')
        except KeyError:
            pass
        trans = []
        mods = dist_info['models']

        for m in mods:
            name = m.pop('model_name')
            mclass = getattr(models, name)
            for p in m:
                if astu.isiterable(m[p]) and 'convert2array' in m[p]:
                    a = np.array(m[p]['value'])
                    a.shape = m[p]['convert2array']
                    m[p] = a
            trans.append(mclass(**m))
        if len(trans) > 1:
            return SCompositeModel(trans, inmap, outmap)
        else:
            return trans[0]

    def undistort(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        if self.distortion is not None:
            return self.distortion(x, y)
        else:
            return x, y

    def __call__(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        if self.distortion is not None:
            x, y = self.distortion(x,y)
        x = self.offx(x)
        y = self.offy(y)
        x, y = self.pcrot(x, y)
        x, y = self.projection(x,y)
        x, y = self.n2c(x, y)
        return x, y

class SpectralWCS(object):
    def __init__(self, spec_info):
        self.spec_info = spec_info
        self.spectral_transform = self.create_spectral_transform()

    def create_spectral_transform(self):
        spec_wcs = copy.deepcopy(self.spec_info)
        model_name = spec_wcs.pop('model_name')
        mclass = getattr(models, model_name)
        return mclass(**spec_wcs)

    def __call__(self, x):
        x = np.asarray(x)
        return self.spectral_transform(x)

class CompositeSpectralWCS(object):
    def __init__(self, wcs_info, spec_info, dist_info=None):
        self.spatial_transform = ImagingWCS(wcs_info, dist_info)
        self.spectral_transform = SpectralWCS(spec_info)

    def undistort(self, x, y):
        input_values = np.asarray(x), np.asarray(y)
        return self.spatial_transform.undistort(x, y)

    def __call__(self, x, y):
        alpha, delta = self.spatial_transform(x, y)
        wave = self.spectral_transform(x)
        return alpha, delta, wave

def get_projcode(ctype):
    return [ctype[0][5:8], ctype[1][5:8]]

def create_projection_transform(projcode):

    projklassname = 'Pix2Sky_' + projcode
    projklass = getattr(projections, projklassname)
    projparams={}
    return projklass(**projparams)

class SelectorModel(Model):
    """
    Parameters
    ----------
    labels : a list of strings or objects
    transforms : a list of transforms
        transforms match labels
    """
    def __init__(self, labels, transforms, n_inputs, n_outputs, param_dim=1):
        self._selector = dict(zip(labels, transforms))
        super(SelectorModel, self).__init__(param_names=[], n_inputs=n_inputs,
                                            n_outputs=n_outputs, param_dim=param_dim)
    @property
    def labels(self):
        return self._selector.keys()

    @property
    def transforms(self):
        returnself._selector.values()

    def __call__(self, label):
        raise NotImplementedError

class RegionsSelector(SelectorModel):
    """
    Parameters
    ----------
    labels : list of strings or ints
        " " and 0 indicate a pixel on the detector which is not within any region
    which is n
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
    def __init__(self, mask_shape, regions_info, wcs_info, wcs_regions=None, spec_regions=None, dist_info=None):
        #create_regions_mask_from_json is defined in the regions API.
        self.regions_mask = region.create_region_mask_from_json(mask_shape, regions_info, validate=False)
        self.reference_region_wcs = wcs_info
        # get the list of region labels
        labels = np.unique(self.regions_mask).tolist()
        try:
            labels.remove(0)
        except ValueError:
            pass
        try:
            labels.remove('')
        except ValueError:
            pass
        transforms = []
        # read in primary WCS and update wcs_relative with it.?
        # read in dist_info and create scomp
        # assign it to rid
        for rid in labels:
            if wcs_regions is not None:
                for reg in wcs_regions:
                    if reg['id'] == rid:
                        rid_wcs = reg
                wcs_info = self.reconstruct_wcs(rid_wcs)
            if dist_info is not None:
                if isinstance(dist_info, list):
                    for reg in dist_info:
                        if reg['id'] == rid:
                            rid_dist_info = reg #['models'][0]
                else:
                    rid_dist_info = dist_info #['models']
            else:
                rid_dist_info = None
            for reg in spec_regions:
                if reg['id'] == rid:
                    spec_wcs = reg['models'][0]
            transform = CompositeSpectralWCS(wcs_info, spec_wcs, rid_dist_info)
            transforms.append(transform)
        super(RegionsSelector, self).__init__(labels, transforms, n_inputs=2, n_outputs=3)

    def reconstruct_wcs(self, rid_wcs):
        wcs_info = copy.deepcopy(self.reference_region_wcs)
        wcs_info['CRVAL'][0] += rid_wcs['CRVAL_DELTA'][0]
        wcs_info['CRVAL'][1] += rid_wcs['CRVAL_DELTA'][1]
        wcs_info['CRPIX'][0] += rid_wcs['CRPIX_DELTA'][0]
        wcs_info['CRPIX'][1] += rid_wcs['CRPIX_DELTA'][1]
        return wcs_info

    def undistort(self, x, y):
        input_values = np.asarray(x), np.asarray(y)
        output_values = [np.zeros_like(arg) for arg in input_values]
        output_values.append(np.zeros_like(input_values[0]))
        result = self.regions_mask[x, y]
        if not astu.isiterable(result):
            if result != 0:
                return self._selector[result[0]](*input_values)
            else:
                return None #input_values
        unique_regions = np.unique(result).tolist()
        unique_regions.remove(0)
        focx = []
        focy = []
        for i in unique_regions:
            indices = (result==i)
            transform = self._selector[i]
            xyind = [val[indices] for val in input_values]
            r, d  = transform.undistort(x[indices], y[indices])
            focx.extend(r)
            focy.extend(d)
        return np.asarray(focx), np.asarray(focy)

    def __call__(self, x, y):
        input_values = np.asarray(x), np.asarray(y)
        output_values = [np.zeros_like(arg) for arg in input_values]
        output_values.append(np.zeros_like(input_values[0]))
        result = self.regions_mask[x, y]
        if not astu.isiterable(result):
            if result != 0:
                return self._selector[result[0]](*input_values)
            else:
                return None #input_values
        unique_regions = np.unique(result).tolist()
        try:
            unique_regions.remove(0)
        except ValueError:
            pass
        try:
            unique_regions.remove("")
        except ValueError:
            pass
        ra = []
        dec = []
        lam = []
        for i in unique_regions:
            indices = (result==i)
            transform = self._selector[i]
            xyind = [val[indices] for val in input_values]
            r, d, l = transform(x[indices], y[indices])
            ra.extend(r)
            dec.extend(d)
            lam.extend(l)
        return np.asarray(ra), np.asarray(dec), np.asarray(lam)

    def inverse(self, x, y):
        """
        Need to invert the regions and construct an InverseSelectorModel
        """

