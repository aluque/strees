from numpy import *
from numpy.linalg import norm
from scipy.special import sph_harm, gammaln

X, Y, Z = 0, 1, 2

def sph_harm_unnorm(m, n, alpha, beta):
    return sph_harm(m, n, alpha, beta) * sqrt(4 * pi / (2 * n + 1))


def expand(q, r, p):
    """ Calculates a multipole expansion or order p centered at the origin
    of a set of charges at points r[3, 0..k-1]. """

    # Using the python convention 0 <= n < p.
    # To access the element n, m we use M[n, m + n].  Hence since the highest
    # m is p - 1, the size of M has to be M[p, 2 * (p - 1) + 1].
    M = zeros((p, 2 * p - 1), dtype='complex128')

    rho = sqrt(sum(r**2, axis=0))
    alpha = arctan2(r[Y, :], r[X, :])
    beta = arccos(r[Z, :] / rho)
    rhon = ones(rho.shape)
    
    for n in xrange(p):
        m = arange(2 * n + 1) - n

        # Y has dimension Y[2n + 1, k]
        Ynm = sph_harm_unnorm(m[:, newaxis], n, alpha, beta).conj()
        M[n, 0: 2 * n + 1] = sum(Ynm * q * rhon, axis=1)
        if n < p - 1:
            rhon = rho * rhon

    return M


def evaluate(r, M):
    phi = zeros((r.shape[1],), dtype='complex128')
    p = M.shape[0]

    rho = sqrt(sum(r**2, axis=0))
    alpha = arctan2(r[Y, :], r[X, :])
    beta = arccos(r[Z, :] / rho)
    rhon = rho
    
    for n in xrange(p):
        m = arange(2 * n + 1) - n
        
        # Y has dimension Y[2n + 1, k]
        Ynm = sph_harm_unnorm(m[:, newaxis], n, alpha, beta)
        T = Ynm.T * M[n, 0: 2 * n + 1] / rhon[:, newaxis]
        
        phi += sum(T, axis=1)
        if n < p - 1:
            rhon = rho * rhon

    return phi


def translate(M, a):
    """ Translates an expansion M centered around a into an expansion
    centered at the origin. """

    p = M.shape[0]

    rho = sqrt(sum(a**2))
    alpha = arctan2(a[Y], r[X])
    beta = arccos(r[Z] / rho)
    
    n = arange(p)

    # Calculate first the factors that depend only on m and n
    for n in arange(p):
        m = arange(2 * n + 1) - n
        F = sph_harm_unnorm(m[:, newaxis], n, alpha, beta).conj()

    for k in xrange(p):
        j = arange(2 * n + 1) - n
        
    
    Ynm = sph_harm_unnorm

def build_A(max_n, max_m):
    n = arange(max_n)
    m = arange(max_m)[newaxis, :]

    return exp(0.5 * (gammaln(n - m + 1) + gammaln(n + m + 1)))

    
def brute(q, r, rdest):
    """ Calculate the potential by 'brute force' i.e. with a k*l algorithm
    Shapes: q[k], r[3, k], rdest[3, l]
    """

    dr = sqrt(sum((r[:, :, newaxis] - rdest[:, newaxis, :])**2, axis=0))
    T = sum(q[:, newaxis] / dr, axis=0)

    return T
    

    
def test():
    import pylab
    
    print build_A(10, 10)
    
    k = 1
    ap = arange(7, 12)
    err = zeros(ap.shape)
    
    r = random.uniform(-1, 1, size=(3, k))
    q = random.uniform(-1.0, 1.0, size=k)
    q[:] = 1.0
    
    X, Y, Z = mgrid[3:5:40j, -1:1:40j, 1::1j]
    x, y = r_[3:5:40j], r_[-1:1:40j]

    rp = c_[X.flat, Y.flat, Z.flat].T
    print rp.shape
    
    #rp = (random.uniform(-0.5, 0.5, size=(3, k * 2))
    #      + array([0.0, 5.0, 5.0])[:, newaxis])

    
    phi_exact = brute(q, r, rp)
    #pylab.plot(w, phi_exact, lw=2.5, c='k')
    pylab.pcolor(x, y, phi_exact.reshape((40, 40)).T)
    pylab.plot(r[0, :], r[1, :], 'o')
    
    for i, p in enumerate(ap):
        M = expand(q, r, p)
        phi = real(evaluate(rp, M))
        err[i] = norm(phi - phi_exact)
        
        print p, err[i]

    pylab.figure(2)
    pylab.pcolor(x, y, phi.reshape((40, 40)).T)
    pylab.plot(r[0, :], r[1, :], 'o')

    pylab.figure(3)
    pylab.plot(ap, err, 'o')
    pylab.show()
    



if __name__ == '__main__':
    test()

