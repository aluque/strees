""" Extract the shape of a branch.
Mostly useful in single-branch simulations."""

from numpy import *
from numpy.linalg import norm
from scipy.optimize import curve_fit
import h5py
import pylab

import datafile


def extract_branch(tr, r, **kwargs):
    """ an array with the distances between segments in the two branches
    after the first branching in a tree. """
    
    ibranch = tr.branches()
    if len(ibranch) > 1:
        raise ValueError("I get confused with more than a single branch")

    ib = ibranch[0]
    
    segment = tr.segments[ib]

    branches = [segment.children[0], segment.children[1]]
    z0 = segment.get(r)[2]
    d = []
    z = []
    
    while True:
        rb = [b.get(r) for b in branches]
        
        if len(branches[0].children) != 1 or len(branches[1].children) != 1:
            break

        branches = [b.children[0] for b in branches]
        z.append(z0 - rb[0][2]) 
        d.append(norm(rb[0] - rb[1]))


    return array(z), array(d)


def analysis(z, y):
    """ Perform some analysis.  Whatever I think of at the moment. """
    print "-" * 20
    # First plot the dat
    pylab.plot(z, y, lw=1.7, c='k')
    pylab.grid(ls='-', c='#999999')

    # asymptote
    dy = y[-1] - y[-100]
    dz = z[-1] - z[-100]

    a = dy / dz
    b = y[-100] - z[-100] * dy / dz

    print "Asymtote y = a * z + b"
    print " a = %f" % a
    print " b = %f" % b

    pylab.plot(z, a*z + b)

    # Fit to an hyperbola
    def f(x, a, b):
        return a * sqrt(x * abs((x - b)))
    
    popt, pcov = curve_fit(f, z, y, [dy / dz, -1.0], maxfev=50000)
    print "Hyperbola"
    print popt

    pylab.plot(z, f(z, *popt))
    

    pylab.show()
    

    
def main():
    """ We called as a stand-alone program, we calculate the angles of the given
    file, step. """
    import sys
    from contextlib import closing

    from optparse import OptionParser
    tfile = None
    params_str = ['external_field', 'branching_sigma', 'tip_mobility']
    parser = OptionParser()

    parser.add_option("--ofile", "-o", dest="ofile", type="str",
                      help="Output file", default=None)

    parser.add_option("--analysis", "-a", dest="analysis", action="store_true",
                      help="Perform some analysis of the data", default=False)

    (opts, args) = parser.parse_args()

    step = args[0]
    files = args[1:]

    for file in files:
        with closing(h5py.File(file)) as fp:
            tr, r = datafile.load_tree(fp, step)
            params = [fp['main'].attrs[k] for k in params_str]
            
        z, d = extract_branch(tr, r)
            
        if opts.ofile is not None:
            savetxt(opts.ofile, c_[z, d])
        
        if opts.analysis:
            analysis(z, d)
            

if __name__ == '__main__':
    main()
    
