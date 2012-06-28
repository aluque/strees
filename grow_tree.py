import sys
import time

from numpy import *
from numpy.random import rand, randn, seed
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
    if RANDOM_SEED >= 0:
        seed(RANDOM_SEED)
    
    EXTERNAL_FIELD_VECTOR = array([0.0, 0.0, EXTERNAL_FIELD])
    
    # init a tree from scratch
    tr, r0, q0 = init_from_scratch()
    
    dt = TIME_STEP
    t = r_[0:END_TIME:dt]
    
    r, q = r0, q0

    dfile = DataFile(OUT_FILE, parameters=parameters)
    branched = False
    
    for i, it in enumerate(t):
        # with ContextTimer("plotting"):
        #     plot_projections(r, q)
        #     pylab.savefig('tree_%.3d.png' % i)
        # print 't = %g\ttree_%.3d.png' % (it, i)
        print "%d/%d  t = %g" % (i, len(t), it)
        branch_prob = BRANCHING_PROBABILITY

        if SINGLE_BRANCHING_TIME > 0:
            if it > SINGLE_BRANCHING_TIME:
                if not branched:
                    branch_prob = inf
                    branched = True

        r, q = step(tr, r, q, it, dt, p=branch_prob)

        with ContextTimer("saving %d" % i):
            dfile.add_step(it, tr, r, q, latest_phi,
                           error=error, error_dq=error_dq)
            

def init_from_scratch():
    """ Init a 'tree' with a single charge point. """
    

    tr = tree.Tree()
    root = tr.make_root()
    r0 = tr.zeros(dim=3)

    # tr = tree.random_branching_tree(4, 0.00)
    # r0 = tree.sample_endpoints(tr) / 1000.0

    k = r0.shape[0]
    q0 = zeros((k, ))

    return tr, r0, q0

    
def step(tr, r, q0, t, dt, p=0.0):
    iterm = tr.terminals()
    box = containing_box(r.T, reflect=HAS_PLANE_ELECTRODE)
    box.set_charges(r.T, q0,
                    max_charges=MAX_CHARGES_PER_BOX,
                    min_length=16 * CONDUCTOR_THICKNESS)
    box.build_lists(recurse=True)
    box.set_field_evaluation(r.T[:, iterm])

    # 1. Calculate the velocities at t
    v0 = velocities(box, tr, r, q0)

    # 2. Relax the tree from t to t + dt
    q1 = relax(box, tr, r, q0, t, dt)

    # 3. Calculate the velocities again at t + dt
    v1 = velocities(box, tr, r, q1)

    # 4. Extend the tree with the leap-frog algo.
    v = 0.5 * (v0 + v1)

    # 5. Branch some of the tips
    vabs = sqrt(sum(v**2, axis=0))
    does_branch = rand(*iterm.shape) < (p * vabs * dt)

    ## if t >= 1.30e-7:
    ##     print "<relaxing only>"
    ##     return r, q1

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


def velocities(box, tr, r, q):
    """ Calculates the electric fields at the tips of the tree and from
    them obtains the propagation velocities of the "streamers" """

    # When we have a single charge the velocity is simply given by the
    # external electri field
    if len(q) == 1:
        return TIP_MOBILITY * EXTERNAL_FIELD_VECTOR[:, newaxis]
    
    box.update_charges(q)
    box.upward(MULTIPOLAR_TERMS)
    box.downward()
    box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
    box.collect_solutions(field=True)
    sfields = self_fields(tr, r, q).T

    v = TIP_MOBILITY * (MAXWELL_FACTOR * box.field
                        + MAXWELL_FACTOR * sfields
                        + EXTERNAL_FIELD_VECTOR[:, newaxis])
    return v

def self_fields(tr, r, q):
    """ Calculates the fields created by the charges at the streamer tips
    on themselves. """
    
    iterm = tr.terminals()
    parents = tr.parents()[iterm]
    
    dr = r[iterm, :] - r[parents, :]
    u = dr / (sqrt(sum(dr**2, axis=1)))[:, newaxis]
    
    return q[iterm][:, newaxis] * u / CONDUCTOR_THICKNESS**2


def relax(box, tr, r, q0, t, dt):
    """ Relax the conductor tree for a time dt. """
    global latest_phi, error, error_dq
    
    #with ContextTimer("re-computing Ohm matrix"):
    # If we have an electrode, we fix q[0] by setting the first row of
    # M to zero.  
    fix = [] if not HAS_PLANE_ELECTRODE else [0]
    
    M = 2 * CONDUCTANCE * tr.ohm_matrix(r, fix=fix)
    n = len(q0)
    
    def f(t0, q):
        global latest_phi, error, error_dq

        if n >= FMM_THRESHOLD:
            # with ContextTimer("FMM") as ct_fmm: 
            box.update_charges(q)

            box.upward(MULTIPOLAR_TERMS)
            box.downward()
            box.solve_all(a=CONDUCTOR_THICKNESS, field=False)

            box.collect_solutions(field=False)

            phi = MAXWELL_FACTOR * box.phi
        else:
            # with ContextTimer("direct") as ct_direct:
            if not HAS_PLANE_ELECTRODE:
                rx, qx = r, q
            else:
                u = array([1, 1, -1])
                rx = concatenate((r, r * u), axis=0)
                qx = r_[q, -q]

            phi0 = MAXWELL_FACTOR * mpolar.direct(rx.T, qx, r.T,
                                                  CONDUCTOR_THICKNESS)
            phi = phi0
            
        # err = sqrt(sum((phi - box.phi)**2)) / len(phi)        
        latest_phi = phi
        error = phi - latest_phi
        error_dq = M.dot(error)
        
        dq = M.dot(phi - dot(r, EXTERNAL_FIELD_VECTOR))

        return dq

    d = ode(f).set_integrator('vode',  nsteps=50000, rtol=1e-7)
    d.set_initial_value(q0, 0.0)
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

