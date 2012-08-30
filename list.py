import sys
import os, os.path
from optparse import OptionParser

from numpy import *
import h5py



def main():
    parser = OptionParser()
    (opts, args) = parser.parse_args()

    for fname in args:
        fp = h5py.File(fname, "r")
        main = fp['main']
        run_name = main.attrs['run_name']
        steps = main.keys()
        
        print "%s [%s]" % (fname, run_name)
        for key, item in main.attrs.iteritems():
            print "\t%-30s =\t%s" % (key, repr(item))

        print "\tLast step: %s [t=%g]" % (steps[-1],
                                          main[steps[-1]].attrs['t'])

    
if __name__ == '__main__':
    main()

    
