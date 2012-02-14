import timeit

from numpy import *
from numpy.linalg import norm
from scipy.special import sph_harm, gammaln

import mpolar
import fmm_3d

X, Y, Z = 0, 1, 2
INWARD, INOUT, OUTWARD = -1, 0, 1
    
def test():
    import pylab
    
    k = 200
    p = 20
    
    r = random.uniform(-1, 1, size=(3, k))
    q = random.uniform(-1.0, 1.0, size=k)

    r0 = array([5.0, 5.0, 5.0])
    r_new = (r + r0[:, newaxis])
    
    reval = array([[1.0, 1.5, 0.5]]).T
    
    M = mpolar.expand(p, r, q, OUTWARD)

    M_new = mpolar.shift((r0[X], r0[Y], r0[Z]), INOUT, M)
    M_direct = mpolar.expand(p, r_new, q, INWARD)

    phi = mpolar.eval_array(M_new, reval, INWARD) 
    phi2 = fmm_3d.brute(q, r_new, reval)
    phi3 = mpolar.direct(r_new, q, reval)
    
    print phi
    print phi2
    print phi3

    print M_new[4, :]
    print M_direct[4, :]
    


if __name__ == '__main__':
    test()

