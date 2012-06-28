from optparse import OptionParser
from glob import iglob
from itertools import cycle, izip

from numpy import *
import pylab
import h5py

X, Y, Z = 0, 1, 2

def main():
    parser = OptionParser()

    (opts, args) = parser.parse_args()

    step = args[0]
    patterns = args[1:]

    colors = ['#888888', '#ff4444', '#4444bb', 'g']
    
    for pattern, color in izip(patterns, cycle(colors)):
        print pattern
        for fname in iglob(pattern):
            print "   ", fname
            fp = h5py.File(fname, "r")
            main = fp['main']
            r = array(main[step]['r'])

            pylab.plot(r[:, X], r[:, Z], 'o', ms=2.0, mec=color, mfc=color),

    pylab.show()

if __name__ == '__main__':
    main()

