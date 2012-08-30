from optparse import OptionParser
from itertools import cycle
from warnings import warn

from numpy import *
import pylab
import h5py

X, Y, Z = 0, 1, 2

def main():
    parser = OptionParser()

    parser.add_option("--parameter", "-p", dest="parameter", action="store",
                      help="Use this parameter", default='conductance')

    parser.add_option("--step", "-s", dest="step", action="store",
                      help="Step to plot", default='last')

    (opts, args) = parser.parse_args()
    fnames = args
    icolors = cycle(['#888888', '#ff4444', '#4444bb', 'g', '#ffbb44',
                     '#ff77ff'])
    d = {}
    
    for fname in fnames:
        fp = h5py.File(fname, "r")
        main = fp['main']
        parm = main.attrs[opts.parameter]
        try:
            steps = main.keys()
        except TypeError:
            warn("Unable to read %s" % fname)
            continue
        
        istep = opts.step if opts.step != 'last' else steps[-2]
        try:
            color = d[parm]
            label = None
        except KeyError:
            color = icolors.next()
            d[parm] = color
            label = "%g" % parm
            
        r = array(main[istep]['r'])
        
        pylab.plot(r[:, X], r[:, Z], 'o', ms=2.0, mec=color, mfc=color,
                   label=label)

    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    main()

