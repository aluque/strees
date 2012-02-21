import sys

from numpy import *
from scipy.integrate import odeint, ode
import pylab

import tree
from refinement import Box, containing_box
import mpolar

from contexttimer import ContextTimer

MAX_CHARGES_PER_BOX = 200
MULTIPOLAR_TERMS = 4
CONDUCTOR_THICKNESS = 0.01
EXTERNAL_FIELD = array([0.0, 0.0, -5.0])
TIP_MOBILITY = 1.0
END_TIME = 1.0
TIME_STEP = 1e-4


def main():
    tr = tree.random_branching_tree(50, 0.05)
    r0 = tree.sample_endpoints(tr) / 1000.0

    k = r0.shape[0]
    q0 = zeros((k, ))

    dt = TIME_STEP
    t = r_[0:END_TIME:dt]
    
    r, q = r0, q0

    for i, it in enumerate(t):
        with ContextTimer("plotting"):
            plot_projections(r, q)
            pylab.savefig('tree_%.3d.png' % i)
        print 't = %g\ttree_%.3d.png' % (it, i)
        
        r, q = step(tr, r, q, dt)
    
    
def step(tr, r, q0, dt):
    iterm = tr.terminals()
    box = containing_box(r.T)
    box.set_charges(r.T, q0, max_charges=MAX_CHARGES_PER_BOX)
    box.build_lists(recurse=True)
    box.set_field_evaluation(r.T[:, iterm])

    # 1. Calculate the velocities at t
    v0 = velocities(box, r, q0)

    # 2. Relax the tree from t to t + dt
    q1 = relax(box, tr, r, q0, dt)

    # 3. Calculate the velocities again at t + dt
    v1 = velocities(box, r, q1)

    # 4. Extend the tree with the leap-frog algo.
    v = 0.5 * (v0 + v1)
    
    radv = r.T[:, iterm] + dt * v
    
    rnew = concatenate((r, radv.T), axis=0)
    qnew = concatenate((q1, zeros((len(iterm), ))), axis=0)
    
    tr.extend(iterm)

    return rnew, qnew


def velocities(box, r, q):
    """ Calculates the electric fields at the tips of the tree and from
    them obtains the propagation velocities of the "streamers" """

    box.update_charges(q)
    box.upward(MULTIPOLAR_TERMS)
    box.downward()
    box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
    box.collect_solutions(field=True)

    return TIP_MOBILITY * (box.field + EXTERNAL_FIELD[:, newaxis])

    
def relax(box, tr, r, q0, dt):
    """ Relax the conductor tree for a time dt. """

    with ContextTimer("re-computing Ohm matrix"):
        M = tr.ohm_matrix(r)

    def f(t0, q):
        with ContextTimer("FMM"): 
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Don't forget to remove these lines:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # box = containing_box(r.T)
            # box.set_charges(r.T, q, max_charges=MAX_CHARGES_PER_BOX)
            # box.build_lists(recurse=True)

            box.clear(recurse=True)
            box.update_charges(q)
            
            box.upward(MULTIPOLAR_TERMS)
            box.downward()
            box.solve_all(a=CONDUCTOR_THICKNESS, field=False)

            box.collect_solutions(field=False)

        # phi = box.phi
        with ContextTimer("direct"):
            phi = mpolar.direct(r.T, q, r.T, CONDUCTOR_THICKNESS)

        err = sum((phi - box.phi)**2) / len(phi)

        print "FMM error = %g" % err
        if err > 1e-4:
            print "max(q) = ", max(q)
            savetxt("poisson.txt", c_[r, q, box.phi, phi])
            sys.exit(-1)
        
        return M.dot(phi - dot(r, EXTERNAL_FIELD))


    d = ode(f).set_integrator('dopri853',  nsteps=2000)
    
    d.set_initial_value(q0, 0.0)

    with ContextTimer("relaxing (n=%d)" % len(q0)):
        d.integrate(dt)

    return d.y


def plot_projections(r, q):
    X, Y, Z = 0, 1, 2
    names = ["X", "Y", "Z"]
    
    axes = [(X, Z), (Y, Z), (X, Y)]
    pylab.subplots_adjust(wspace=0.35, hspace=0.25, right=0.95, top=0.95)
    for i, (d1, d2) in enumerate(axes):
        ax = pylab.subplot(2, 2, i)

        ax.clear()
        ax.scatter(r[:, d1], r[:, d2], c=q,
                      s=5.0, faceted=False),

        ax.set_xlabel(names[d1])
        ax.set_ylabel(names[d2])
        
        #ax.quiver(r[iterm, 0], r[iterm, 2], field[0, :], field[2, :])
        #ax.set_xlim([0.0, 1.0])
        #ax.set_ylim([-0.2, 1.0])


if __name__ == '__main__':
    main()

