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

    t = zeros((len(steps)),)
    r = zeros((len(steps)),)

    try:
        for i, step in enumerate(steps):
            ri = array(main[step]['r'])

            t[i] = main[step].attrs['t']
            r[i] = amax(sqrt(sum(ri**2, axis=1)))

            print "%g\t%g\t#%s" % (t[i], r[i], step)

    except KeyError:
        pass

    dt = t[1] - t[0]
    v = diff(r) / dt

    tmid = 0.5 * (t[1:] + t[:-1])

    #a, b = simple_regression(log(tmid[-100:]), log(v[-100:]))
    #print a
    
    if opts.show:
        pylab.figure(1)
        
        pylab.xlabel("t")
        pylab.ylabel("r")
        pylab.plot(t, r, lw=1.7)

        pylab.figure(2)
        pylab.xlabel("t")
        pylab.ylabel("v")
        pylab.plot(tmid, v, lw=1.7)
        #pylab.plot(tmid, exp(b + a * log(tmid)))
        pylab.loglog()
        pylab.show()
        

    if opts.ofile is not None:
        savetxt(opts.ofile, c_[t, r])


def simple_regression(xi, yi):
    A = ones((len(yi), 2), dtype=float)
    A[:, 0] = xi[:]
    
    r = linalg.lstsq(A, yi)
    return r[0][0], r[0][1]



    
if __name__ == '__main__':
    main()

    
