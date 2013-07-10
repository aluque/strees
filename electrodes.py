""" Implementation of electrode geometries with the image charge method.
"""
from numpy import *

import tree
X, Y, Z = 0, 1, 2

class Electrode(object):
    def __init__(self):
        pass

    def images(self, r, q):
        """ Calculates the images of charges q located at points r.
        Returns rimag, qimag. """

        raise NotImplementedError

    def extend(self, r, q):
        rimag, qimag = self.images(r, q)
        
        return concatenate((r, rimag), axis=0), r_[q, qimag]
    
    
class NullElectrode(object):
    """ Simply a trick to simulate no electrode. """
    def images(self, r, q):
        return zeros((0, 3)), zeros((0,))
    

class Planar(Electrode):
    def __init__(self, x0, axis=Z):
        """ Inits a plane electrode perpendicular to the axis axis, located
        at x0. """
        self.x0 = x0
        self.u = ones((3, ))
        self.u[axis] *= -1

    def images(self, r, q):
        return self.images_r(r), self.images_q(r, q)
    
    def images_q(self, r, q):
        return -q

    def images_r(self, r):
        return 2 * self.x0 + r * self.u


class Sphere(Electrode):
    def __init__(self, center, a):
        """ A spherical electrode centered at center and with radius a. """
        self.center = center
        self.a = a
        self.a2 = a * 2
        
    def images(self, r, q):
        # See e.g. Jackson, 2.2
        dr = r - self.center[newaxis, :]
        y = sqrt(sum(dr**2, axis=1))
        
        yp = self.a2 / y
        qp = -dist.q * self.a / y
        rp = self.center[newaxis, :] + dr * yp / y

        return rp, qp


    def images_q(self, r, q):
        dr = r - self.center[newaxis, :]
        y = sqrt(sum(dr**2, axis=1))
        
        return (-q * self.a / y)

        
    def images_r(self, r):
        dr = r - self.center[newaxis, :]
        y = sqrt(sum(dr**2, axis=1))
        
        yp = self.a2 / y
        rp = self.center[newaxis, :] + dr * yp / y

        return rp
