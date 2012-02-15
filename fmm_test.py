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

    # r = array([[0.3, 0.2, 0.1]]).T
    
    # q = array([1.0])
    
    reval = array([[2.0, 3.0, 1.0]]).T
    
    M = mpolar.expand(p, r, q, OUTWARD)
    print M
    
    ef = mpolar.eval_field_array(M, reval, OUTWARD) 
    ef_0 = mpolar.field_direct(r, q, reval, 0.0)
    
    print ef, norm(ef)
    print ef_0, norm(ef_0)
    
    print M[:, 0]

if __name__ == '__main__':
    test()

