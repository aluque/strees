import hotshot, hotshot.stats
import time

from numpy import *
import pylab
from matplotlib import cm

import tree
import mpolar
from refinement import Box


def main():
    p = 5
    k = 50000

    t0 = time.time()
    
    tr = tree.random_branching_tree(k, 0.05)
    r = tree.sample_endpoints(tr)
    
    r0 = amin(r, axis=0)
    r1 = amax(r, axis=0)
    
    r = (r - r0) / (r1 - r0)
    r0 = amin(r, axis=0)
    r1 = amax(r, axis=0)
    
    q = ones((r.shape[0],))
    
    t1 = time.time()

    # Let's play with trees
    parents = tr.parents()
    l = sqrt(sum((r - r[parents, :])**2, axis=1))
    print l.shape

    box = Box(r0, r1)
    box.set_charges(r.T, q, max_charges=200)

    box.build_lists(recurse=True)
    box.upward(p)
    box.downward()
    box.solve_all()
    box.collect_solutions()

    t2 = time.time()

    phi = mpolar.direct(r.T, q, r.T, 0.0)

    t3 = time.time()

    savetxt("cmp.txt", c_[phi, box.phi])
    
    eps = sqrt(sum((phi - box.phi)**2) / (k * (k - 1)))
    
    print ("t(setup) = %g s  t(multipol) = %g s  t(direct) = %g s"
           % (t1 - t0, t2 - t1, t3 - t2))
    print "TOTAL: %g" % (t3 - t0)
    print "ERROR: %g" % eps

    
    # box.plot(dims=[0, 2], recurse=True, fill=False)

    # plot_tree(tr, r, box.phi)
    # pylab.figure(2)
    # plot_tree(tr, r, phi)
    # pylab.show()


def plot_tree(tr, r, phi):
    vmin, vmax = amin(phi), amax(phi)
    cmap = cm.get_cmap('jet')
    #cmap.set_clim(vmin=vmin, vmax=vmax)

    for segment in tr:
        ep = segment.get(r)
        try:
            ip = segment.parent.get(r)
        except AttributeError:
            ip = ep
        c = cmap((segment.get(phi) - vmin) / (vmax - vmin))
        pylab.plot([ip[0], ep[0]], [ip[2], ep[2]], lw=1.2, c=c)




if __name__ == '__main__':
    prof = hotshot.Profile("solve_tree.prof")
    prof.runcall(main)
    prof.close()
    stats = hotshot.stats.load("solve_tree.prof")
    stats.strip_dirs()
    stats.sort_stats('time', 'calls')
    stats.print_stats(20)

