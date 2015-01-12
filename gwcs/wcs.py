# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import numpy as np
from astropy.io import fits
from astropy.modeling import models
from astropy.utils import isiterable

from . import coordinate_systems
from .util import ModelDimensionalityError, CoordinateFrameError
from .selector import *


__all__ = ['WCS']


class WCS(object):
    """
    Basic WCS class

    Parameters
    ----------
    output_coordinate_system : str, gwcs.coordinate_systems.Frame
        A coordinates object or a string label
    input_coordinate_system : str, gwcs.coordinate_systems.Frame
        A coordinates object or a string label
    forward_transform : astropy.modeling.Model
        a model to do the forward transform
    name : str
        a name for this WCS
    """
    def __init__(self, output_coordinate_system,  input_coordinate_system='detector',
                 forward_transform=None, name=""):
        self._forward_transform = forward_transform
        self._input_coordinate_system = input_coordinate_system
        self._output_coordinate_system = output_coordinate_system
        self._name = name
        '''
        if forward_transform is not None and input_coordinate_system is not None \
           and output_coordinate_system is not None:
            self._pipeline.add_transform(self._input_coordinate_system.__class__,
                                         self._output_coordinate_system.__class__,
                                         forward_transform)
        '''

    @property
    def unit(self):
        return self._output_coordinate_system._unit

    @property
    def output_coordinate_system(self):
        return self._output_coordinate_system

    @property
    def input_coordinate_system(self):
        return self._input_coordinate_system

    @property
    def forward_transform(self):
        return self._forward_transform

    @forward_transform.setter
    def forward_transform(self, value):
        self._forward_transform = value.copy()

    @property
    def as_world_coordinates(self, *args):
        return self.output_coordinate_system.world_coordinates(*args)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    def __call__(self, *args):
        if self._forward_transform is not None:
            result = self._forward_transform(*args)
        if self.output_coordinate_system is not None:
            return self.output_coordinate_system.world_coordinates(*result)
        else:
            return result

    #def invert(self, *args, **kwargs):
        #'''
        #args : x, y
        #**kwargs: x0, y0, maxiter ,...
        #'''
        #try:
            #return self.forward_transform.inverse(*args)
        #except (NotImplementedError, KeyError):
            ##if initial_values == {}:
            #if len(args) == 2:
                #try:
                    #foc2world = self.forward_transform['focal_plane' : ]
                    #x0, y0 = foc2world.inverse(*args)
                    #args += (x0, y0)
                #except IndexError:
                    #raise CoordinateFrameError("Coordinate frame 'focal_plane' not found."
                                               #"Either provide initial values or assign a "
                                               #"label 'focal_plane'")
                #except NotImplementedError:
                    #import warnings
                    #warnings.warn("Cannot compute initial values, Using defaults (0, 0).")
                    #initial_values = {'x0': 0.0, 'y0': 0.0}
                    #x0, y0 = (initial_values['x0'], initial_values['y0'])
                    #args += (x0, y0)
            #else:
                ##if 'x0' not in initial_values or 'y0' not in initial_values:
                #if len(args) != 4:
                    #raise ValueError("Expected `initial values` to be a dictionary "
                                     #"with keys `x0` and `y0`.")
            #return self._invert(*args, **kwargs)

    #def _invert(self, x, y, x0, y0, **kwargs):
        #maxiter = kwargs.pop('maxiter', 20)
        #accuracy = kwargs.pop('accuracy', 1.0e-4)
        #adaptive = kwargs.pop('adaptive', False)
        #detect_divergence = kwargs.pop('detect_divergence', True)
        #quiet = kwargs.pop('quiet', False)
        #if isiterable(x0):
            #x = x0.copy()
            #y = y0.copy()
        #else:
            #x = x0
            #y = y0

        #foc2world = self.forward_transform['focal_plane' : ]
        #dx, dy = foc2world(x, y)
        #dx -= x0
        #dy -= y0

        ## update initial solution:
        #x -= dx
        #y -= dy

        ## norm (L2) squared of the correction:
        #dn2prev = dx**2+dy**2
        #dn2 = dn2prev
        ## prepare for iterative process
        #iterlist = range(1, maxiter+1)
        #print('iterlist', iterlist)
        #accuracy2 = accuracy**2
        #ind = None
        #inddiv = None
        #old_invalid = np.geterr()['invalid']
        #old_over = np.geterr()['over']
        #np.seterr(invalid = 'ignore', over = 'ignore')

        ######################################################################
        ###                     NON-ADAPTIVE ITERATIONS:                    ##
        ######################################################################
        #if not adaptive:
            #for k in iterlist:
                ## check convergence:
                #if np.max(dn2) < accuracy2:
                    #break

                ## find correction to the previous solution:
                #'''
                #dx, dy = self.pix2foc(x, y, origin)
                #'''
                #dx, dy = foc2world(x, y)
                ## If pix2foc does not apply all the required distortion
                ## corrections then replace the above line with:
                ##r0, d0 = self.all_pix2world(x, y, origin)
                ##dx, dy = self.wcs_world2pix(r0, d0, origin )
                #dx -= x0
                #dy -= y0

                ## update norn (L2) squared of the correction:
                #dn2 = dx**2+dy**2

                ## check for divergence (we do this in two stages
                ## to optimize performance for the most common
                ## scenario when succesive approximations converge):
                #if detect_divergence:
                    #print('dn2, dn2prev', dn2, dn2prev)
                    #ind, = np.where(dn2 <= dn2prev)
                    #if ind.shape[0] < npts:
                        #inddiv, = np.where(
                            #np.logical_and(dn2 > dn2prev, dn2 >= accuracy2))
                        #if inddiv.shape[0] > 0:
                            ## apply correction only to the converging points:
                            #if np.isscalar(x):
                                #x -= dx
                                #y -= dy
                            #else:
                                #x[ind] -= dx[ind]
                                #y[ind] -= dy[ind]

                            ## switch to adaptive iterations:
                            #ind, = np.where((dn2 >= accuracy2) & \
                                            #(dn2 <= dn2prev) & np.isfinite(dn2))
                            #iterlist = iterlist[k:]
                            #adaptive = True
                            #break
                    ##dn2prev[ind] = dn2[ind]
                    #dn2prev = dn2

                ## apply correction:
                #x -= dx
                #y -= dy
         ######################################################################
        ###                      ADAPTIVE ITERATIONS:                       ##
        ######################################################################
        #if adaptive:
            #if ind is None:
                #ind = np.asarray(range(npts), dtype=np.int64)

            #for k in iterlist:
                ## check convergence:
                #if ind.shape[0] == 0:
                    #break

                ## find correction to the previous solution:
                #'''
                #dx[ind], dy[ind] = self.pix2foc(x[ind], y[ind], origin)
                #'''
                #if isiterable(dx):
                    #dx[ind], dy[ind] = foc2world( x[ind], y[ind])

                    ## If pix2foc does not apply all the required distortion
                    ## corrections then replace the above line with:
                    ##r0[ind], d0[ind] = self.all_pix2world(x[ind], y[ind], origin)
                    ##dx[ind], dy[ind] = self.wcs_world2pix(r0[ind], d0[ind], origin)
                    #dx[ind] -= x0[ind]
                    #dy[ind] -= y0[ind]
                #else:
                    #dx, dy = foc2world(x, y)
                    #dx -=x0
                    #dy -= y0
                ## update norn (L2) squared of the correction:
                #dn2 = dx**2+dy**2

                ## update indices of elements that still need correction:
                #if detect_divergence:
                    #ind, = np.where((dn2 >= accuracy2) & (dn2 <= dn2prev))
                    ##ind = ind[np.where((dn2[ind] >= accuracy2) & (dn2[ind] <= dn2prev))]
                    #dn2prev[ind] = dn2[ind]
                #else:
                    #ind, = np.where(dn2 >= accuracy2)
                    ##ind = ind[np.where(dn2[ind] >= accuracy2)]

                ## apply correction:
                #if isiterable(dx):
                    #x[ind] -= dx[ind]
                    #y[ind] -= dy[ind]
                #else:
                    #x -= dx
                    #y -= dy
         ######################################################################
        ###         FINAL DETECTION OF INVALID, DIVERGING,                  ##
        ###         AND FAILED-TO-CONVERGE POINTS                           ##
        ######################################################################
        ## Identify diverging and/or invalid points:
        #invalid = (((~np.isfinite(y)) | (~np.isfinite(x)) | \
                    #(~np.isfinite(dn2))) & \
                   #(np.isfinite(ra)) &  (np.isfinite(dec)))
        ## When detect_divergence==False, dn2prev is outdated (it is the
        ## norm^2 of the very first correction). Still better than nothing...
        #inddiv, = np.where(((dn2 >= accuracy2) & (dn2 > dn2prev)) | invalid)
        #if inddiv.shape[0] == 0:
            #inddiv = None
        ## identify points that did not converge within
        ## 'maxiter' iterations:
        #if k >= maxiter:
            #ind,= np.where((dn2 >= accuracy2) & (dn2 <= dn2prev) & (~invalid))
            #if ind.shape[0] == 0:
                #ind = None
        #else:
            #ind = None
        ######################################################################
        ###      RAISE EXCEPTION IF DIVERGING OR TOO SLOWLY CONVERGING      ##
        ###      DATA POINTS HAVE BEEN DETECTED:                            ##
        ######################################################################
        ## raise exception if diverging or too slowly converging
        #if (ind is not None or inddiv is not None) and not quiet:
            #print('ind', ind, inddiv)
            ##if vect1D:
            #sol  = [x, y]
            #err  = [np.abs(dx), np.abs(dy)]
            #'''
            #else:
                #sol  = np.dstack( [x, y] )[0]
                #err  = np.dstack( [np.abs(dx), np.abs(dy)] )[0]
            #'''
            ## restore previous numpy error settings:
            #np.seterr(invalid = old_invalid, over = old_over)

            #if inddiv is None:
                #raise NoConvergenceError("'WCS._invert()' failed to "
                                    #"converge to the requested accuracy after {:d} "         \
                                    #"iterations.".format(k), best_solution = sol,            \
                                    #accuracy = err, niter = k, failed2converge = ind,        \
                                    #divergent = None)
            #else:
                #raise NoConvergenceError("'WCS._invert()' failed to "
                                    #"converge to the requested accuracy.{0:s}"               \
                                    #"After {1:d} iterations, the solution is diverging "     \
                                    #"at least for one input point."                          \
                                    #.format(os.linesep, k), best_solution = sol,             \
                                    #accuracy = err, niter = k, failed2converge = ind,        \
                                    #divergent = inddiv)
         ######################################################################
        ###             FINALIZE AND FORMAT DATA FOR RETURN:                ##
        ######################################################################
        ## restore previous numpy error settings:
        #np.seterr(invalid = old_invalid, over = old_over)

        ##if vect1D:
        #return x, y
        ###else:
        ###    return np.dstack( [x, y] )[0]

    def transform(self, fromsys, tosys, *args):
        """
        Perform coordinate transformation beteen two frames inclusive.

        Parameters
        ----------
        fromsys : CoordinateFrame
            an instance of CoordinateFrame
        tosys : CoordinateFrame
            an instance of CoordinateFrame
        args : float
            input coordinates to transform
        """
        transform = self._forward_transform[fromsys : tosys]
        return transform(*args)


    def get_transform(self, fromsys, tosys):
        """
        Return a transform between two coordinate frames

        Parameters
        ----------
        fromsys : CoordinateFrame
            an instance of CoordinateFrame
        tosys : CoordinateFrame
            an instance of CoordinateFrame
        """
        try:
            return self._forward_transform[fromsys : tosys]
        except ValueError:
            try:
                transform = self._forward_transform[tosys : fromsys]
            except ValueError:
                return None
            try:
                return transform.inverse
            except NotImplementedError:
                return None

    @property
    def available_frames(self):
        return self.forward_transform.submodel_names

    def footprint(self, axes, center=True):
        """
        Parameters
        ----------
        axes : tuple of floats
            size of image
        center : bool
            If `True` use the center of the pixel, otherwise use the corner.

        Returns
        -------
        coord : (4, 2) array of (*x*, *y*) coordinates.
            The order is counter-clockwise starting with the bottom left corner.
        """
        naxis1, naxis2 = axes # extend this to more than 2 axes
        if center == True:
            corners = np.array([[1, 1],
                                [1, naxis2],
                                [naxis1, naxis2],
                                [naxis1, 1]], dtype = np.float64)
        else:
            corners = np.array([[0.5, 0.5],
                                [0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, 0.5]], dtype = np.float64)
        return self.__call__(corners[:,0], corners[:,1])
        #result = np.vstack(self.__call__(corners[:,0], corners[:,1])).T
        #try:
            #return self.output_coordinate_system.world_coordinates(result[:,0], result[:,1])
        #except:
            #return result

