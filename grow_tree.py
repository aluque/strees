""" grow_tree.py is the main module of the code and the one you invoke from
the command line to start a simulation. 

To start a simulation with parameters read from a file simulation.ini
simply invoke this module as::

   python grow_tree.py simulation.ini

The output will be written in a file name simulation.h5 with a 
`HDF5 <http://www.hdfgroup.org/HDF5//>`_ format.

"""
# NOTES:
#
# In this module and in all the rest, the geometrical dimension always
# corresponds to the last axis of an array.  For example to save the locations
# of k points we use an array with shape [k, 3].
#
# All-caps variable names are in general global variables that are read
# from the input file.  We use globals().update(...) so they are avaliable
# everywhere.
#

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
import electrodes
from contexttimer import ContextTimer

latest_phi = None
EXTERNAL_FIELD_VECTOR = None
ELECTRODE = None
X, Y, Z = 0, 1, 2

class TooLongTimestep(Exception):
    pass

def main():
    """ This is the main function of the code it is the starting point of
    a simulation. """

    # Load input parameters from the input file and add the, in allcaps
    # to the global namespace.
    global EXTERNAL_FIELD_VECTOR, ELECTRODE
    parameters = load_input(sys.argv[1], param_descriptors)
    globals().update(dict((key.upper(), item)
                          for key, item in parameters.iteritems()))
    if RANDOM_SEED >= 0:
        seed(RANDOM_SEED)
    
    EXTERNAL_FIELD_VECTOR = array([0.0, 0.0, EXTERNAL_FIELD])
    ELECTRODE = init_electrode()
    
    # init a tree from scratch
    tr, r0, q0 = init_from_scratch(INITIAL_NODES)
    
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

        if SINGLE_BRANCHING_Z != 0 and not branched:
            zterm = r[tr.terminals()[0], Z]
            if zterm < SINGLE_BRANCHING_Z:
                if not branched:
                    branch_prob = inf
                    branched = True

        r, q = adapt_step(tr, r, q, dt, p=branch_prob)

        with ContextTimer("saving %d" % i):
            phi = solve_phi(r, q)
            dfile.add_step(it, tr, r, q, phi,
                           error=error, error_dq=error_dq)
            
        if END_WITH_RECONNECTION and tr.reconnects(r):
            print "Finishing due to a reconnection."
            break
    

def init_from_scratch(n=0):
    """ Init a 'tree' with the root node plus n additional nodes in a vertical 
    string. """
    
    tr = tree.Tree()
    root = tr.make_root()

    for i in xrange(n):
        tr.extend([i,])

    r0 = tr.zeros(dim=3)

    k = r0.shape[0]

    r0[:, Z] = -arange(k) * CONDUCTOR_THICKNESS

    q0 = zeros((k, ))

    return tr, r0, q0

    
def adapt_step(tr, r0, q0, dt, p=0.0):
    """ Performs a step of duration dt but divides it into sub steps
    to make sure that the length of a channel is never longer than ``MAX_STEP``.
    """
    current_dt = dt
    r, q = r0, q0
    remaining_steps = 1
    while remaining_steps > 0:
        try:
            r, q = step(tr, r, q, current_dt, p=p) 
            remaining_steps -= 1
        except TooLongTimestep:
            current_dt /= 2.
            remaining_steps *= 2
    return r, q


def step(tr, r, q0, dt, p=0.0):
    """ Performs an elementary step, including relaxation and advancing
    the channels. 

    Arguments:

      * *tr*: the :class:`tree.Tree` instance containing the tree structure.
      * *r*: an array containing the node locations.
      * *q0*: an array containing the charges of the nodes.
      * *dt*: the time step.

    """

    iterm = tr.terminals()

    box = containing_box(r, electrode=ELECTRODE)
    box.set_charges(r, q0,
                    max_charges=MAX_CHARGES_PER_BOX,
                    min_length=16 * CONDUCTOR_THICKNESS)
    box.build_lists(recurse=True)
    box.set_field_evaluation(r[iterm, :])

    # 1. Calculate the velocities at t
    v0 = velocities(box, tr, r, q0)

    # 2. Relax the tree from t to t + dt
    q1 = relax(box, tr, r, q0, dt)

    # 3. Calculate the velocities again at t + dt
    v1 = velocities(box, tr, r, q1)

    # 4. Extend the tree with the leap-frog algo.
    v = 0.5 * (v0 + v1)
    
    # 5. Branch some of the tips
    vabs = sqrt(sum(v**2, axis=1))

    # If the longest step is longer than MAX_STEP, raise an exception
    # telling the calling function to reduce dt.
    if (max(vabs) * dt) > MAX_STEP:
        raise TooLongTimestep

    does_branch = rand(*iterm.shape) < (p * vabs * dt)

    radv = empty((sum(does_branch) + sum(vabs > 0), 3))
    j = 0
    
    for i, branches in enumerate(does_branch):
        if not branches:
            if vabs[i] > 0:
                radv[j, :] = r[iterm[i], :] + dt * v[i, :]
                j += 1
        else:
            # Note that slow channels, although unlikely, may branch.
            # However, not if their velocity is 0
            dr1, dr2 = symmetric_gaussian(dt * v[i, :], BRANCHING_SIGMA)
            radv[j, :] = r[iterm[i], :] + dr1
            radv[j + 1, :] = r[iterm[i], :] + dr2
            j += 2
    
    rnew = concatenate((r, radv), axis=0)
    qnew = concatenate((q1, zeros((sum(does_branch) 
                                   + sum(vabs > 0),))), axis=0)
    
    tr.extend(sort(r_[iterm[vabs > 0], 
                      iterm[does_branch]]))
    return rnew, qnew


def velocities(box, tr, r, q):
    """ Calculates the electric fields at the tips of the tree and from
    them obtains the propagation velocities of the *streamers* """

    iterm = tr.terminals()

    # When we have a single charge the velocity is simply given by the
    # external electric field
    if len(q) == 1:
        return TIP_MOBILITY * external_field(r[iterm, :])
    
    box.update_charges(q)
    box.upward(MULTIPOLAR_TERMS)
    box.downward()
    box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
    box.collect_solutions(field=True)
    sfields = self_fields(tr, r, q)
    E = (MAXWELL_FACTOR * box.field
         + MAXWELL_FACTOR * sfields
         + external_field(r[iterm, :]))

    absE = sqrt(sum(E**2, axis=1))

    # An unit vector with the same direction as E
    u = E / absE[:, newaxis]

    # Now we can calculate the absolute value of the velocity
    vabs = TIP_MOBILITY * where(absE > TIP_MIN_FIELD, absE - TIP_MIN_FIELD, 0)

    v = u * vabs[:, newaxis]
    
    return v


def self_fields(tr, r, q):
    """ Calculates the fields created by the charges at the streamer tips
    on themselves. """
    
    iterm = tr.terminals()
    parents = tr.parents()[iterm]
    
    dr = r[iterm, :] - r[parents, :]
    u = dr / (sqrt(sum(dr**2, axis=1)))[:, newaxis]
    
    return q[iterm][:, newaxis] * u / CONDUCTOR_THICKNESS**2


def relax(box, tr, r, q0, dt):
    """ Relax the conductor :class:`tree.Tree` *tr* for a time *dt*. 

    Arguments:

      * *tr*: the :class:`tree.Tree` instance containing the tree structure.
      * *r*: an array containing the node locations.
      * *q0*: an array containing the charges of the nodes.
      * *dt*: the time step.
    """
    global latest_phi, error, error_dq
    
    #with ContextTimer("re-computing Ohm matrix"):
    # If we have an electrode, we fix q[0] by setting the first row of
    # M to zero.  
    fix = [] if ELECTRODE is None else [0]
    
    # On Fri Aug 31 11:46:47 2012 I found a factor 2 here that I do not know
    # where it comes from.  Probably it was a reliq of the mid-points approach
    # (But it was duplicated in ohm_matrix anyway!).  I am removing it.
    M = CONDUCTANCE * tr.ohm_matrix(r, fix=fix)
    n = len(q0)
    
    def f(t0, q):
        global latest_phi, error, error_dq

        phi = solve_phi(r, q, box)
        # err = sqrt(sum((phi - box.phi)**2)) / len(phi)        
        latest_phi = phi
        error = phi - latest_phi
        error_dq = M.dot(error)
        
        dq = M.dot(phi + external_potential(r))

        return dq

    d = ode(f).set_integrator('vode',  nsteps=250000, rtol=1e-8)
    d.set_initial_value(q0, 0.0)
    d.integrate(dt)

    return d.y


def solve_phi(r, q, box=None):
    if len(q) >= FMM_THRESHOLD and box is not None:
        # with ContextTimer("FMM") as ct_fmm: 
        box.update_charges(q)

        box.upward(MULTIPOLAR_TERMS)
        box.downward()
        box.solve_all(a=CONDUCTOR_THICKNESS, field=False)

        box.collect_solutions(field=False)

        phi = MAXWELL_FACTOR * box.phi
    else:
        # with ContextTimer("direct") as ct_direct:
        if ELECTRODE is None:
            rx, qx = r, q
        else:
            rx, qx = ELECTRODE.extend(r, q)

        phi0 = MAXWELL_FACTOR * mpolar.direct(rx, qx, r,
                                              CONDUCTOR_THICKNESS)
        phi = phi0

    return phi


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

    if not BRANCH_IN_XZ:
        p, q = sigma * randn(2)
    else:
        if FIXED_BRANCHING_ANGLE > 0:
            p, q = norm(dr) * tan(FIXED_BRANCHING_ANGLE / 2), 0.0
        else:
            p, q = sigma * randn(), 0.0
    
    dr1 = dr + (p * e1 + q * e2)
    dr2 = dr - (p * e1 + q * e2)

    if FIXED_BRANCHING_ANGLE > 0:
        # This is to avoid too long segments at branching points.
        # Presently I am doing it only here to preserve compatibility
        # with the algorithm described in the paper as of Sat Mar 23 20:58:46 2013
        dr1 *= norm(dr) / norm(dr1)
        dr2 *= norm(dr) / norm(dr2)

    return dr1, dr2


def external_field(r):
    """ Calculates the external field at points *r*.  This is calculated
    from ``EXTERNAL_FIELD`` and ``ELECTRODE_POTENTIAL``.  As the code stands now
    only these two possibilities are physically meaningful:

     1. Specify ``EXTERNAL_FIELD`` with a planar electrode or with no electrode,
        but use ``ELECTRODE_POTENTIAL=0``.
     2. ``ELECTRODE_POTENTIAL != 0``, but ``ELECTRODE_GEOMETRY = 'sphere'`` and
        ``EXTERNAL_FIELD = 0``.

    However, we allow the user to shoot himself on his foot, so he can
    select any arbitrary combination of these parameters.  Beware.
    """
    field = EXTERNAL_FIELD_VECTOR[newaxis, :]
    
    if ELECTRODE_POTENTIAL == 0:
        return field

    center = array([0.0, 0.0, ELECTRODE_RADIUS])
    dr = r - center[newaxis, :]
    rabs = sqrt(sum(dr**2, axis=1))
    field = field + (ELECTRODE_RADIUS * ELECTRODE_POTENTIAL
                     * dr / rabs[:, newaxis]**3)

    return field


def external_potential(r):
    """ Calculates the external potential at points *r*.  See above, in
    external_field for the risks here.
    """
    phi = -dot(r, EXTERNAL_FIELD_VECTOR)

    if ELECTRODE_POTENTIAL == 0:
        return phi

    center = array([0.0, 0.0, ELECTRODE_RADIUS])
    dr = r - center[newaxis, :]
    rabs = sqrt(sum(dr**2, axis=1))
    phi = phi + ELECTRODE_RADIUS * ELECTRODE_POTENTIAL / rabs

    return phi



def init_electrode():
    """ Uses the input parameters to select an electrode geometry. """
    def planar():
        """ Planar electrode.  Always located at z=0. """
        return electrodes.Planar(0)

    def sphere():
        """ Sphere electrode.  Located at [0, 0, -ELECTRODE_RADIUS]. """
        center = array([0.0, 0.0, ELECTRODE_RADIUS])
        return electrodes.Sphere(center, ELECTRODE_RADIUS)

    def null():
        """ No electrode. """
        # return electrodes.NullElectrode()
        # This is actually faster:
        return None

    d = dict(planar=planar, plane=planar, sphere=sphere, null=null, none=null)
    
    try:
        return d[ELECTRODE_GEOMETRY]()
    except KeyError:
        raise KeyError("Electrode geometry '%s' not recognized"
                       % ELECTRODE_GEOMETRY)


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

