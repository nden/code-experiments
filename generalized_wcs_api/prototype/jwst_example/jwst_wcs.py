# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import copy
import numpy as np
from astropy.modeling import *
from astropy import utils as astu
from . import region
import json

class JwstWcs(WCS):
    """
    WCS for all JWST instruments.

    Note: The data is in fits format but reference files are in json.

    Parameters
    ----------
    meta : json
        a storage for meta data for an observation, basicaaly a dict of {keyword: value}
    header : astropy.fits.Header
        fits header
    dist_info : json
        reference file with distortion models
    spec_wcs : json
        reference file with spectral models
    regions_def : json
        reference file with definitions of IFU regions as polygon vertices in detector coordinates
    wcs_regions : json
        reference file with basic WCS for each IFU region
    spec_regions : json
        reference file with spectral models for each IFU region
    """

    def __init__(self, meta, header, dist_info=None, spec_wcs=None, regions_def=None,
                 wcs_regions=None, spec_regions=None):
        self.mode = get_observing_mode(header)
        self.wcs_info = self.read_wcs_from_header(header, self.mode)
        self.distortion = dist_info
        self.regions_def = regions_def
        self.wcs_regions = wcs_regions
        self.spec_regions = spec_regions
        dateobs = meta.observation.date.isoformat()
        coord_system = self.create_coordinate_system(self.wcs_info, dateobs)
        if self.mode == "imaging":
            transform = transforms.ImagingTransform(self.wcs_info, dist_info)#wcs_transform(self.wcs_info, dist_info)
        elif self.mode == 'spectroscopic':
            if len(self.wcs_info['CTYPE']) == 3:
                transform = transforms.CompositeSpectralTransform(self.wcs_info, spec_wcs['models'][0], dist_info)
            else:
                print(self.wcs_info['CTYPE'])
                transform = transforms.SpectralTransform(spec_wcs['models'][0])
        elif self.mode == 'ifu':
            assert regions_def is not None,"regions definition file is required for IFU observations"
            assert wcs_regions is not None, "wcs_regions reference file is required for IFU observations"
            assert spec_regions is not None, "spec_regions reference file is required for IFU observations"
            mask_shape = (meta.subarray.xsize, meta.subarray.ysize)
            transform = transforms.RegionsSelector(mask_shape, regions_def,
                                self.wcs_info, wcs_regions, spec_regions, dist_info )
        elif self.mode == 'multislit':
            assert regions_def is not None, "regions definition file is required for IFU observations"
            assert spec_regions is not None, "spec_regions reference file is required for IFU observations"
            mask_shape = (meta.subarray.xsize, meta.subarray.ysize)
            transform = transforms.RegionsSelector(mask_shape, regions_def,
                                self.wcs_info, wcs_regions, spec_regions, dist_info )

        super(JwstWcs, self).__init__(transform, output_coordinate_system=coord_system)

    def read_wcs_from_header(self, header, mode):
        """
        Read basic WCS info from the Primary header of the data file.
        """
        wcs_info = {}
        if mode == 'imaging':
            wcsaxes = 2
        elif mode == 'spectroscopic':
            wcsaxes = 3
            #wcs_info['SPECSYS'] = header['SPECSYS']
        elif mode == 'ifu':
            wcsaxes = 3
        elif mode == 'multislit':
            wcsaxes = 3
        else:
            raise ValueError('Unrecognized mode')
        WCSAXES = header.get('WCSAXES', wcsaxes)
        wcs_info['WCSAXES'] = WCSAXES
        wcs_info['RADESYS'] = header.get('RADESYS', 'ICRS')
        wcs_info['VAFACTOR'] = header.get('VAFACTOR', 1)
        wcs_info['CUNIT'] = [header.get('CUNIT'+str(i), "deg") for i in range(1, WCSAXES+1)]#3)]
        wcs_info['CRPIX'] = [header.get('CRPIX'+str(i), 0.) for i in range(1, 3)]
        wcs_info['CRVAL'] = [header.get('CRVAL'+str(i), 0.) for i in range(1, WCSAXES+1)]#3)]
        wcs_info['CTYPE'] = [header.get('CTYPE'+str(i), "") for i in range(1, WCSAXES+1)]#3)]
        wcs_info['CDELT'] = [header.get('CDELT'+str(i), 1.) for i in range(1, WCSAXES+1)]#3)]
        wcs_info['PC'] = [header.get('PC'+str(i)+'_'+str(j)) for i in range(1, 3) for j in range(1, 3)]
        return wcs_info

    def undistort(self, x, y):
        return self.transform.undistort(x, y)
        #if self.transform.distortion is not None:
            #return self.transform.distortion(x, y)
        #else:
            #return x, y

    def get_projcode(self, ctype):
        return [ctype[0][5:8], ctype[1][5:8]]

    def create_coordinate_system(self, wcs_info, dateobs):
        """
        This is almost a stub function.
        In real the headers of spectroscopic observations should have
        WCSAXES = 3
        CTYPE3 = WAVE

        or

        WCSAXES =1 ?

        if the spatial information is not important.
        """
        ref_system = wcs_info['RADESYS']
        unit = wcs_info['CUNIT']
        ctype = wcs_info['CTYPE']
        names = [label.split('-')[0] for label in ctype]
        projcode = self.get_projcode(ctype)[0].upper()
        return coordinate_systems.SkyFrame(ref_system, reference_position="BARYCENTER",
                                           obstime=dateobs,#time.Time(dateobs, scale='TT', format='iso'),
                                           projection=projcode,
                                           axes_names=[ctype[0].split('-')[0],
                                                       ctype[1].split('-')[0]],
                                           units=["deg", "deg"])


class ImagingTransform(object):
    """
    Implements a composite transform: distortion + WCS + projection

    """
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

class SpectralTransform(object):
    
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

class CompositeSpectralTransform(object):
    
    def __init__(self, wcs_info, spec_info, dist_info=None):
        self.spatial_transform = ImagingTransform(wcs_info, dist_info)
        self.spectral_transform = SpectralTransform(spec_info)

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

def get_observing_mode(header):
        """                                                                                      Currently all example data files have a 'MODE'  keyword
        but in real this will probably look for the EXP_TYPE keyword.
        """
        try:
            exp_type = header['EXP_TYPE']
        except KeyError:
            if header['INSTRUME'].lower() == 'nirspec':
                if header['SUBARRAY'].lower() in ['allslits', 'full']:
                    mode = 'MULTISLIT'
                else:
                    raise

        if exp_type in ['MIR_LRS','NIRISS_SOSS']:
            mode = 'spectroscopic'
        elif exp_type in ['MIR_MRS','NRS_IFU']:
            mode = 'ifu'
        elif exp_type in ['NRS_FIXEDSLIT','NRS_MSA']:
            mode = 'multislit' 
        else:
            mode = 'imaging'

        return mode.lower()
