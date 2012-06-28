""" Implementation of the charge simulation technique. """

import sys

from numpy import *
import pylab
import h5py
from scipy.sparse.linalg import LinearOperator, bicgstab, bicg, gmres

from refinement import Box, containing_box
import mpolar
from plotter import bounding_box, plot_projections

def paraboloid(a, ztip, nz, ntheta, phase_shift=False):
    """ Generates points in a paraboloid z - z0 = a * r**2 with ntheta points
    in each circular cross-section located at points z.  If phase_shift is
    True the points are phase-shifted pi / ntheta radians. """

    theta = linspace(0, 2 * pi, ntheta)[:-1]
    z = linspace(ztip, 0, nz + 1)[1:-1]

    if phase_shift:
        theta = theta + (theta[1] - theta[0]) / 2
    
    rho = sqrt((z - ztip) / a)

    rx = rho * cos(theta[:, newaxis])
    ry = rho * sin(theta[:, newaxis])
    rz = z + 0 * theta[:, newaxis]

    return c_[r_[0, ravel(rx)], r_[0, ravel(ry)], r_[ztip, ravel(rz)]].T



def direct_with_electrode(r, q, reval, conductor_thickness):
    u = array([1, 1, -1])[:, newaxis]
    rx = concatenate((r, r * u), axis=1)
    qx = r_[q, -q]

    phi = mpolar.direct(rx, qx, reval, conductor_thickness)

    return phi


def build_linop(r, reval, conductor_thickness):
    def f_phi(q):
        return direct_with_electrode(r, q, reval, conductor_thickness)

    return LinearOperator((reval.shape[1], r.shape[1]), f_phi, dtype='float64')


def main():
    ztip = -0.01
    
    rp = paraboloid(800.0, ztip, 16, 16)
    rp_phi = paraboloid(1.0, ztip, 16, 16, phase_shift=True)
    np = rp.shape[1]

    print rp.shape, rp_phi.shape
    fname, step = sys.argv[1:3]

    fp = h5py.File(fname, "r")
    g = fp['main']
    r = array(g[step]['r']).T
    q = array(g[step]['q'])
    phi0 = array(g[step]['phi'])

    CONDUCTOR_THICKNESS = fp['main'].attrs['conductor_thickness']
    MAX_CHARGES_PER_BOX = fp['main'].attrs['max_charges_per_box']
    MULTIPOLAR_TERMS = fp['main'].attrs['multipolar_terms']
    HAS_PLANAR_ELECTRODE = fp['main'].attrs['has_plane_electrode']
    print "%d charges" % len(q)

    r[2, :] += ztip
    rt = concatenate((rp, r), axis=1)

    # 1. Compute the potential created at the sources and in the phi points
    #    of the electrode
    phi = direct_with_electrode(r, q, rt, CONDUCTOR_THICKNESS)
    phip = phi[:np]
    print phip
    
    # 2. Solve the charges that will make phi zero in the paraboloid points.
    linop = build_linop(rp, rp_phi, CONDUCTOR_THICKNESS)
    qp, info = gmres(linop, -phip, tol=1e-6)
    print qp
    
    r0, r1 = bounding_box(rt)
    plot_projections(rt.T, r_[qp, 0 * q], r0, r1, plot3d=True)
    pylab.show()



if __name__ == '__main__':
    main()
    
