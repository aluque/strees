import sys
import os, os.path
from optparse import OptionParser
from matplotlib.colors import LogNorm

from numpy import *
import h5py

try:
    import pylab
except ImportError:
    pass

X, Y, Z = 0, 1, 2

def main():
    parser = OptionParser()

    parser.add_option("--show", dest="show", action="store_true",
                      help="Open the matplotlib window?", default=False)

    parser.add_option("--ofile", "-o", dest="ofile", action="store",
                      help="Save to this file", default=None)

    (opts, args) = parser.parse_args()

    fname = args[0]

    fp = h5py.File(fname, "r")
    main = fp['main']
    run_name = main.attrs['run_name']
    steps = main.keys()

    eta = main.attrs['conductance']
    mu_T = main.attrs['tip_mobility']
    E0 = main.attrs['external_field']
    external_field_vector = array([0.0, 0.0, E0])

    eps0 = 1 / (4 * pi * main.attrs['maxwell_factor'])
    d = main.attrs['conductor_thickness']

    q0 = eta**2 / (mu_T**2 * E0 * eps0)
    ell = eta / (mu_T * E0 * eps0)
    tau = ell / (mu_T * E0)
    
    t = zeros((len(steps)),)
    q = zeros((len(steps)),)
    r = zeros((len(steps)),)
    E = zeros((len(steps)),)
    
    for i, step in enumerate(steps):
        qi = array(main[step]['q'])
        ri = array(main[step]['r'])
        phi = array(main[step]['phi'])
        phi = phi - dot(ri, external_field_vector)

        t[i] = main[step].attrs['t']
        q[i] = sum(qi)
        r[i] = weighted_r(ri, qi)
        E[i] = (phi[1] - phi[0]) / abs(ri[0, Z] - ri[1, Z])

        print "%g\t%g\t%g\t%g\t#%s" % (t[i], q[i], r[i], E[i], step)

    if opts.show:
        pylab.xlabel("t")
        pylab.ylabel("Q")
        pylab.plot(t, q, lw=1.7)
        pylab.show()

    if opts.ofile is not None:
        savetxt(opts.ofile, c_[t, q, r, E, t / tau, q / q0, r / ell])


def weighted_r(r, q):
    """ Finds the average |r| weighted with q. """
    rmod = sqrt(sum(r**2, axis=1))
    try:
        return average(rmod, weights=q)
    except ZeroDivisionError:
        # If there are no charges...
        return 0.0


if __name__ == '__main__':
    main()

    
