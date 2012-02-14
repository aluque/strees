import time

from numpy import *
from scipy.sparse import lil_matrix, csr_matrix
from scipy.integrate import odeint, ode
import pylab
from matplotlib import cm

import tree
import mpolar
from refinement import Box


def main():
    p = 5
    k = 5000

    tr = tree.random_branching_tree(k, 0.05)
    r = tree.sample_endpoints(tr)
    
    r0 = amin(r, axis=0)
    r1 = amax(r, axis=0)
    
    r = (r - r0) / (r1 - r0)
    r0 = amin(r, axis=0)
    r1 = amax(r, axis=0)
    
    # Let's play with trees
    parents = tr.parents()
    l = sqrt(sum((r - r[parents, :])**2, axis=1))
    rmid = 0.5 * (r + r[parents, :])
    
    M = tree_matrix(tr, l)

    t = linspace(0, 1e-1, 1000)
    q0 = zeros((r.shape[0],))
    f = build_func(M, rmid.T, r0, r1,
                   array([0.0, 0.0, -1.0]),
                   a=0.01, for_ode=True)
    
    # qt = odeint(f, q0, t, atol=1e-1)
    d = ode(f).set_integrator('dopri5', method='bdf')

    d.set_initial_value(q0, 0.0)
    dt = t[1] - t[0]
    qt = zeros((len(t) + 1, len(q0)))
    i = 0
    while d.successful() and d.t < t[-1]:
        i += 1
        d.integrate(d.t + dt)
        print d.t
        qt[i, :] = d.y
        
    animate(tr, rmid, qt)


def build_func(M, r, r0, r1, e0, a=0.0, p=4, for_ode=False):
    def f(q, t0):
        box = Box(r0, r1)
        box.set_charges(r, q, max_charges=200)

        box.build_lists(recurse=True)
        box.upward(p)
        box.downward()
        box.solve_all(a=a)

        box.collect_solutions()

        return M.dot(box.phi - dot(r.T, e0))

    def f_ode(t0, q):
        return f(q, t0)

    return f if not for_ode else f_ode


def animate(tr, r, qt):
    nt, n = qt.shape
    vmin, vmax = amin(qt), amax(qt)
    
    for i in xrange(nt):
        pylab.clf()
        #plot_tree(tr, r, qt[i, :], vmin=vmin, vmax=vmax)
        print r[:, 0].shape, r[:, 2].shape, qt[i, :].shape
        
        pylab.scatter(r[:, 0], r[:, 2], c=qt[i, :],
                      s=5.0, faceted=False,
                      vmin=vmin, vmax=vmax)
        pylab.xlim([0.0, 1.0])
        pylab.ylim([0.0, 1.0])
        
        pylab.savefig('tree_%.3d.png' % i)
        print 'tree_%.3d.png' % i


def tree_matrix(tr, l):
    """ Builds a matrix M that will give use the evolution of charges
    in every segment of the tree as dq/dt = M . phi, where phi is
    the potential at the center of each segment and '.' is the dot product.
    """
    n = l.shape[0]

    # We build the matrix in LIL format first, later we convert to a
    # format more efficient for matrix-vector multiplications
    M = lil_matrix((n, n))

    for segment in tr:
        i = segment.index
        m = 0.0
        for other in segment.iter_adjacent():
            j = other.index
            a = 2.0 / (l[i] + l[j])
            M[i, j] = a
            m -= a

        M[i, i] = m
        
    return csr_matrix(M)


def plot_tree(tr, r, phi, vmin=None, vmax=None):
    if vmin is None:
        vmin = amin(phi)

    if vmax is None:
        vmax = amax(phi)
        
    cmap = cm.get_cmap('jet')
    #cmap.set_clim(vmin=vmin, vmax=vmax)

    for segment in tr:
        ep = segment.get(r)
        try:
            ip = segment.parent.get(r)
        except AttributeError:
            ip = ep
        #print segment.get(phi)
        c = cmap((segment.get(phi) - vmin) / (vmax - vmin))
        
        pylab.plot([ip[0], ep[0]], [ip[2], ep[2]], lw=1.2, c=c)


if __name__ == '__main__':
    main()
