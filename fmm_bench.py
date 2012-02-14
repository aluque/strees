import timeit

from numpy import *
from numpy.linalg import norm
from scipy.special import sph_harm, gammaln

import mpolar
import fmm_3d

X, Y, Z = 0, 1, 2

    
def test():
    import pylab
    
    k = 200
    p = 6
    
    r = random.uniform(-1, 1, size=(3, k))
    q = random.uniform(-1.0, 1.0, size=k)

    reval = array([[5.0, 5.0, 5.0]]).T + random.uniform(-1, 1, size=(3, k))
    
    def py_func():
        M = fmm_3d.expand(q, r, p)
        phi = fmm_3d.evaluate(reval, M)
        #print M
        return M
    
    def c_func():
        M = mpolar.expand(p, r, q, 1)
        phi = mpolar.eval_array(M, reval, 1)
        #print M
        return M

    def brute_func():
        phi = fmm_3d.brute(q, r, reval)
    
    def direct_func():
        phi = mpolar.direct(r, q, reval)
        
    M_py = py_func()
    M_c = c_func()
    brute_func()

    # for l in xrange(p):
    #     #print M_py[l, l: 2 * l + 1]
    #     #print M_c[l, :l + 1]
    #     f = M_py[l, l: 2 * l + 1] / M_c[l, :l + 1]
    #     print f**2
        
        
    for f in (py_func, c_func, brute_func, direct_func):
        print timeit.timeit(stmt=f, number=10)
        


if __name__ == '__main__':
    test()

