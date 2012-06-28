""" Calculation of the branching angles of a tree. """

from numpy import *
from numpy.linalg import norm
import h5py

import datafile

def branching_angles(tr, r, **kwargs):
    """ Returns an array with the branching angles of the tree. """
    
    ibranch = tr.branches()
    angles = zeros(ibranch.shape)

    for i, ib in enumerate(ibranch):
        angles[i] = angle1(tr.segments[ib], r, **kwargs)

    return angles


def angle1(segment, r, max_n=15, skip_n=0):
    """ Gets the angle that forms the branch at segment.  Sees further by
    max_n segments.  If there is another branch inside these segments,
    measures the angle only up to that branching point. """

    branches = [segment.children[0], segment.children[1]]
    origin = segment.get(r)
    
    for i in xrange(max_n):
        rb = [b.get(r) for b in branches]
        if i == skip_n:
            r0 = rb
        
        if len(branches[0].children) != 1 or len(branches[1].children) != 1:
            break

        branches = [b.children[0] for b in branches]
        
    u = [ri - r0i for ri, r0i in zip(rb, r0)]

    return arccos(dot(u[0], u[1]) / norm(u[0]) / norm(u[1]))
    

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

    parser.add_option("--degrees", "-d", dest="degrees", action="store_true",
                      help="Use sexadecimal degrees instead of radians",
                      default=False)

    parser.add_option("--max", "-m", dest="max_n", action="store",
                      type="int", 
                      help="Look at N segments below the branching point",
                      default=15)

    parser.add_option("--skip", "-s", dest="skip_n", action="store",
                      type="int", 
                      help="Ignore N segments below the branching point",
                      default=0)

    parser.add_option("--table", "-t", dest="tfile", type="str",
                      help="Write a table in a file", default=None)

    (opts, args) = parser.parse_args()

    step = args[0]
    files = args[1:]

    if opts.tfile is not None:
        tfile = open(opts.tfile, "a")
    
    for file in files:
        with closing(h5py.File(file)) as fp:
            tr, r = datafile.load_tree(fp, step)
            params = [fp['main'].attrs[k] for k in params_str]
            
        a = branching_angles(tr, r,
                             max_n=opts.max_n,
                             skip_n=opts.skip_n)

        if opts.degrees:
            a[:] = a * 180 / pi

        if opts.ofile is not None:
            savetxt(opts.ofile, a)

        avg = mean(a)
        sigma = std(a)

        if tfile:
            tfile.write("\t".join(str(x) for x in (params + [avg, sigma]))
                        + '\t# %s\n' % (file))
            

        print "[%s] mean = %g\tstd deviation = %g" % (file, avg, sigma)

    if tfile:
        tfile.close()

if __name__ == '__main__':
    main()
    
