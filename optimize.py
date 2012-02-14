""" Time a few runs of the adaptative FMM to try to optimize its parameters.
"""

import timeit

from numpy import *
import mpolar
from refinement import Box
import tree

def run_instance(k, max_charges, p, number=1):
    """ Runs an instance of the FMM with k charges, max_charges in each box
    and multipolar order p.  Returns the time required for the run. """

    # For charges in a tree:
    #
    tr = tree.random_branching_tree(k, 0.05)
    r = tree.sample_endpoints(tr)
    
    r0 = amin(r, axis=0)
    r1 = amax(r, axis=0)
    
    r = ((r - r0) / (r1 - r0)).T
    r0 = amin(r, axis=1)
    r1 = amax(r, axis=1)
    
    q = ones((r.shape[1],))
    
    # For uniformly distributed charges:
    #
    # r = random.uniform(-1, 1, size=(3, k))
    # q = random.uniform(-1.0, 1.0, size=k)

    # r0 = array([-1.0, -1.0, -1.0])
    # r1 = array([1.0, 1.0, 1.0])

    def mp_func():
        """ Fast Multipolar Method """
        box = Box(r0, r1)
        box.set_charges(r, q, max_charges=max_charges)
        box.build_lists(recurse=True)
        box.upward(p)
        box.downward()
        box.solve_all()
        box.collect_solutions()


    return timeit.timeit(stmt=mp_func, number=number)

def direct_time(k, number=1):
    r = random.uniform(-1, 1, size=(3, k))
    q = random.uniform(-1.0, 1.0, size=k)

    def direct_func():
        """ Direct method """
        return mpolar.direct(r, q, r)

    return timeit.timeit(stmt=direct_func, number=number)



def main():
    k = 10000
    p = 4
    
    max_charges = exp(linspace(log(10), log(10000), 50)).astype('i')
    
    deltas = zeros(max_charges.shape)
    for i, mc in enumerate(max_charges):
        print "max_charges = %d" % mc
        delta = run_instance(k, mc, p)
        print "Duration = %f" % delta
        deltas[i] = delta

    direct = zeros(deltas.shape) + direct_time(k)
    
    savetxt('opt_k%d_p%d_tree.txt' % (k, p), c_[max_charges, deltas, direct])


if __name__ == '__main__':
    main()
    
        
