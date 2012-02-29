import sys
import time

from numpy import *
from numpy.random import rand, randn
from numpy.linalg import norm
from scipy.integrate import odeint, ode

from readinput import load_input
import tree
from refinement import Box, containing_box
import mpolar
from datafile import DataFile
import os.path
from angles import branching_angles
import parameters as param_descriptors

from contexttimer import ContextTimer

# MAX_CHARGES_PER_BOX = 200
# MULTIPOLAR_TERMS = 5
# CONDUCTOR_THICKNESS = 0.001
# EXTERNAL_FIELD = array([0.0, 0.0, -5.0])
# TIP_MOBILITY = 1.0
# END_TIME = 1.0
# TIME_STEP = 1e-4
# BRANCHING_PROBABILITY = 100.
# BRANCHING_SIGMA = 1e-4
# RUN_NAME = sys.argv[1]
# OUT_FILE = os.path.expanduser('~/data/trees/%s.h5' % sys.argv[1])

# # Use the FMM solver only when we have more than this number of charges
# FMM_THRESHOLD = 4000

latest_phi = None
EXTERNAL_FIELD_VECTOR = None

def main():
    # Load input parameters from the input file and add the, in allcaps
    # to the global namespace.
    global EXTERNAL_FIELD_VECTOR
    parameters = load_input(sys.argv[1], param_descriptors)
    globals().update(dict((key.upper(), item)
                          for key, item in parameters.iteritems()))

    EXTERNAL_FIELD_VECTOR = array([0.0, 0.0, EXTERNAL_FIELD])
    
    tr = tree.random_branching_tree(50, 0.05)
    r0 = tree.sample_endpoints(tr) / 1000.0

    k = r0.shape[0]
    q0 = zeros((k, ))

    dt = TIME_STEP
    t = r_[0:END_TIME:dt]
    
    r, q = r0, q0


    dfile = DataFile(OUT_FILE, parameters=parameters)
    
    for i, it in enumerate(t):
        # with ContextTimer("plotting"):
        #     plot_projections(r, q)
        #     pylab.savefig('tree_%.3d.png' % i)
        # print 't = %g\ttree_%.3d.png' % (it, i)

        r, q = step(tr, r, q, dt, p=BRANCHING_PROBABILITY)
        #with ContextTimer("saving"):
        dfile.add_step(it, tr, r, q, latest_phi)
            

    
def step(tr, r, q0, dt, p=0.0):
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
    
    # 5. Branch some of the tips
    does_branch = rand(*iterm.shape) < (p * dt)

    # print "%d active branches" % len(iterm)
    
    radv = empty((3, sum(does_branch) + len(iterm)))
    j = 0
    for i, branches in enumerate(does_branch):
        if not branches:
            radv[:, j] = r.T[:, iterm[i]] + dt * v[:, i]
            j += 1
        else:
            dr1, dr2 = symmetric_gaussian(dt * v[:, i], BRANCHING_SIGMA)
            radv[:, j] = r.T[:, iterm[i]] + dr1
            radv[:, j + 1] = r.T[:, iterm[i]] + dr2
            j += 2
    
    rnew = concatenate((r, radv.T), axis=0)
    qnew = concatenate((q1, zeros((sum(does_branch) + len(iterm), ))), axis=0)
    
    tr.extend(sort(r_[iterm, iterm[does_branch]]))

    return rnew, qnew


def velocities(box, r, q):
    """ Calculates the electric fields at the tips of the tree and from
    them obtains the propagation velocities of the "streamers" """

    box.update_charges(q)
    box.upward(MULTIPOLAR_TERMS)
    box.downward()
    box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
    box.collect_solutions(field=True)
    
    return TIP_MOBILITY * (box.field + EXTERNAL_FIELD_VECTOR[:, newaxis])

    
def relax(box, tr, r, q0, dt):
    """ Relax the conductor tree for a time dt. """
    global latest_phi
    
    #with ContextTimer("re-computing Ohm matrix"):
    M = tr.ohm_matrix(r)

    n = len(q0)
    
    def f(t0, q):
        global latest_phi

        if n >= FMM_THRESHOLD:
            # with ContextTimer("FMM") as ct_fmm: 
            box.update_charges(q)

            box.upward(MULTIPOLAR_TERMS)
            box.downward()
            box.solve_all(a=CONDUCTOR_THICKNESS, field=False)

            box.collect_solutions(field=False)

            phi = box.phi
        else:
            # with ContextTimer("direct") as ct_direct:
            phi = mpolar.direct(r.T, q, r.T, CONDUCTOR_THICKNESS)

        # err = sqrt(sum((phi - box.phi)**2)) / len(phi)
        
        # ftimes.write("%d\t%g\t%g\t%g\n" % (len(q),
        #                                   ct_fmm.duration, ct_direct.duration,
        #                                   err))
        # ftimes.flush()
        
        # print "FMM error = %g" % err
        latest_phi = phi
        
        return M.dot(phi - dot(r, EXTERNAL_FIELD_VECTOR))


    d = ode(f).set_integrator('dopri853',  nsteps=2000)
    
    d.set_initial_value(q0, 0.0)

    #with ContextTimer("relaxing (n=%d)" % len(q0)):
    d.integrate(dt)

    return d.y


def symmetric_gaussian(dr, sigma):
    """ Samples a branch from a symmetric, gaussian branching model.
    In a plane perpendicular to dr we sample dr1 from a cylindrically
    symmetric gaussian distribution; the two branching points are dr1 and
    its symmetric vector wrt dr. """

    u = dr / norm(dr)
    # We find two unit vectors orthonormal to u (also dr); note that this
    # fails if u is parallel to x !!!
    
    ex = array([1.0, 0, 0])

    e1 = ex - dot(u, ex) * u
    e1 = e1 / norm(e1)
    
    e2 = cross(u, e1)

    p, q = sigma * randn(2)

    dr1 = dr + (p * e1 + q * e2)
    dr2 = dr - (p * e1 + q * e2)

    return dr1, dr2


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

