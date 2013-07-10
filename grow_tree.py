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
    global SPRITE_RS, SPRITE_QS
    parameters = load_input(sys.argv[1], param_descriptors)
    globals().update(dict((key.upper(), item)
                          for key, item in parameters.iteritems()))
    if RANDOM_SEED >= 0:
        seed(RANDOM_SEED)
    
    EXTERNAL_FIELD_VECTOR = array([0.0, 0.0, EXTERNAL_FIELD])
    ELECTRODE = init_electrode()

    # Initialization of sprite fields
    if SPRITES:
        SPRITE_RS, SPRITE_QS = sprite_charge_sources()

    # init a tree from scratch
    tr, dist0 = init_from_scratch(INITIAL_NODES)
    
    dt = TIME_STEP
    t = r_[0:END_TIME:dt]
    
    dist = dist0

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

        dist = adapt_step(tr, dist, it, dt, p=branch_prob)

        with ContextTimer("saving %d" % i):
            # Mon Jul  1 20:41:57 2013: COMPAT. BREAKING.
            # From now on, I will store the full potential, not
            # only the self-constistent potential, which only had
            # sense when debugging.
            phi = solve_phi(dist) + external_potential(dist.r, it)
            dfile.add_step(it, tr, dist, phi,
                           error=error, error_dq=error_dq)
            
        if END_WITH_RECONNECTION and tr.reconnects(dist):
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
    dist = tree.Distribution(r=r0, q=q0)

    return tr, dist

    
def adapt_step(tr, dist0, t, dt, p=0.0):
    """ Performs a step of duration dt but divides it into sub steps
    to make sure that the length of a channel is never longer than ``MAX_STEP``.
    """
    current_dt = dt
    current_t = t
    dist = dist0
    remaining_steps = 1
    while remaining_steps > 0:
        try:
            dist = step(tr, dist, current_t, current_dt, p=p) 
            print "t = %g; dt = %g" % (current_t, current_dt)
            remaining_steps -= 1
            current_t += current_dt
        except TooLongTimestep:
            current_dt /= 2.
            remaining_steps *= 2
    return dist


def step(tr, dist, t, dt, p=0.0):
    """ Performs an elementary step, including relaxation and advancing
    the channels. 

    Arguments:

      * *tr*: the :class:`tree.Tree` instance containing the tree structure.
      * *r*: an array containing the node locations.
      * *q0*: an array containing the charges of the nodes.
      * *t*: Initial time of the time-step.  This is only used for time-
             dependent fields.
      * *dt*: the time step.

    """

    iterm = tr.terminals()

    box = containing_box(dist.r, electrode=ELECTRODE)
    box.set_charges(dist.r, dist.q,
                    max_charges=MAX_CHARGES_PER_BOX,
                    min_length=16 * CONDUCTOR_THICKNESS)
    box.build_lists(recurse=True)
    box.set_field_evaluation(dist.r[iterm, :])

    # 1. Calculate the velocities at t
    v0 = velocities(box, tr, dist, t)

    # 2. Relax the tree from t to t + dt
    dist1 = relax(box, tr, dist, t, dt)

    # 3. Calculate the velocities again at t + dt
    v1 = velocities(box, tr, dist, t + dt)

    # 4. Extend the tree with the leap-frog algo.
    v = 0.5 * (v0 + v1)
    
    # 5. Advance the tips and branch some of them
    dist_new = advance(tr, dist1, v, dt, p=p, iterm=iterm)
    
    # 6. Add emerging upward streamers (for sprites)
    if SPRITES:
        dist_new = upward_streamers(tr, dist_new, t, dt)

    return dist_new



def advance(tr, dist, v, dt, p=0.0, iterm=None):
    """ Advances the tree by adding new nodes, one per each terminal
    node if it is not branching, two if it branches.

    Aguments:

      * *tr*: the :class:`tree.Tree` instance containing the tree structure.
              Note that this will be updated to allow for the newly created
              nodes.
      * *dist*: The distribution of r, q, etc in the tree.
      * *v*: an array containing the velocities of the terminal nodes.
      * *dt*: the time step.
      * *p*: The branching rate.
      * *iterm*: If given, the indices of the terminal nodes. If none, 
                 they will be calculated from *tr*.
    
    """

    if iterm is None:
        iterm = tr.terminals()

    vabs = sqrt(sum(v**2, axis=1))

    if not SPRITES:
        theta = 1
    else:
        theta = sprite_theta(dist.r[iterm, :])

    # If the longest step is longer than MAX_STEP, raise an exception
    # telling the calling function to reduce dt.
    if (max(vabs) * dt > MAX_STEP / theta).any():
        raise TooLongTimestep

    does_branch = rand(*iterm.shape) < (theta * p * vabs * dt)

    radv = empty((sum(does_branch) + sum(vabs > 0), 3))
    j = 0
    
    for i, branches in enumerate(does_branch):
        if not branches:
            if vabs[i] > 0:
                radv[j, :] = dist.r[iterm[i], :] + dt * v[i, :]
                j += 1
        else:
            # Note that slow channels, although unlikely, may branch.
            # However, not if their velocity is 0
            factor = 1 if not SPRITES else 1 / theta[i]
            dr1, dr2 = symmetric_gaussian(dt * v[i, :], 
                                          factor * BRANCHING_SIGMA)
            radv[j, :] = dist.r[iterm[i], :] + dr1
            radv[j + 1, :] = dist.r[iterm[i], :] + dr2
            j += 2
    
    rnew = concatenate((dist.r, radv), axis=0)
    qnew = concatenate((dist.q, zeros((sum(does_branch) 
                                       + sum(vabs > 0),))), axis=0)
    
    tr.extend(sort(r_[iterm[vabs > 0], 
                      iterm[does_branch]]))
    return tree.Distribution(r=rnew, q=qnew)

has_upward = False
def upward_streamers(tr, dist, t, dt):
    """ Decides if new upward-propagating streamers will emerge from the
    body of the streamer tree.  This is a key point for the simulation
    of carrot sprites. """
    
    global has_upward

    if has_upward:
        return dist

    p = tr.parents()
    l = tr.lengths(dist)

    # We consider that the charge contained in one segment is half of
    # each of the charges at its two extremes.
    qsegment = 0.5 * (dist.q + dist.q[p])
    lmbd = qsegment / l

    if not SPRITES:
        theta = 1
    else:
        midpoints = tr.midpoints(dist)
        theta = sprite_theta(midpoints)

    # The radial field created by the charges in the tree.
    er = MAXWELL_FACTOR * lmbd * theta / (2 * CONDUCTOR_THICKNESS)
    e0 = TRANS_BREAK_FIELD * theta
    t0 = TRANS_BREAK_TIME / theta
    l0 = TRANS_BREAK_LENGTH / theta

    # This is the rate of transversal breakdown at each segment
    # Now we allow trans. breakdown only from negatively charged segments.
    # Yes, one can drop theta above, but it helps me seeing them there
    w = l * (abs(er / e0) ** TRANS_BREAK_ALPHA - 1) / t0 / l0
    w[:] = where(er < 0, w, 0)

    print "max(er / e0) = {:g}".format(nanmax(where(er < 0, abs(er/e0), 0)))

    ibranches = nonzero(rand(*w.shape) < w * dt)[0]

    if len(ibranches) == 0:
        # Nothing to do here
        return dist

    radv = empty((len(ibranches), 3))
    field0 = mpolar.field_direct_ap(dist.r, dist.q, dist.r[ibranches, :], 
                                    CONDUCTOR_THICKNESS / theta)
    field = MAXWELL_FACTOR * field0 + external_field(dist.r[ibranches, :], t)

    for i, ib in enumerate(ibranches):
        # Skip the terminal nodes.
        if not tr.segments[ib].children:
            continue

        # This would generate a random perp. initial direction:
        e1, e2 = perp_basis(dist.r[ib] - dist.r[p[ib]])
        phi = 2 * pi * rand()
        enew = e1 * cos(phi) + e2 * sin(phi)
        vabs = TIP_MOBILITY * er[ib]
        if SPRITES:
            vabs /= theta[ib]
        radv[i, :] = dist.r[ib] + enew * vabs * dt

        # This is to follow the (-) local electric field but does not look into
        # the repulsion from the channel
        # radv[i, :] = r[ib] - TIP_MOBILITY * field[i] * dt
        print "New upward branch!"
        print "r[ib] = ", dist.r[ib]
        print "radv[i] = ", radv[i, :] 

        has_upward = True


    rnew = concatenate((r, radv), axis=0)
    qnew = concatenate((q, zeros((len(ibranches),))), axis=0)
    tr.extend(ibranches)

    return tree.Distribution(r=rnew, q=qnew)


def velocities(box, tr, dist, t):
    """ Calculates the electric fields at the tips of the tree and from
    them obtains the propagation velocities of the *streamers* """

    iterm = tr.terminals()

    # When we have a single charge the velocity is simply given by the
    # external electric field
    if len(dist.q) == 1:
        return TIP_MOBILITY * external_field(dist.r[iterm, :], t)
    
    if len(dist.q) >= FMM_THRESHOLD and box is not None:
        box.update_charges(dist.q)
        box.upward(MULTIPOLAR_TERMS)
        box.downward()
        box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
        box.collect_solutions(field=True)
        field = box.field
    else:
        if not SPRITES:
            field = mpolar.field_direct(dist.r, dist.q, dist.r[iterm, :], 
                                        CONDUCTOR_THICKNESS)
        else:
            theta = sprite_theta(dist.r)
            field = mpolar.field_direct_ap(dist.r, dist.q, dist.r[iterm, :], 
                                           CONDUCTOR_THICKNESS / theta)

    sfields = self_fields(tr, dist)
    # The sign(q) appears because we allow negative streamers to propagate
    # upwards.
    E = (sign(dist.q[iterm][:, newaxis]) 
         * (MAXWELL_FACTOR * field + external_field(dist.r[iterm, :], t))
         + MAXWELL_FACTOR * sfields)

    flt = dist.q[iterm] < 0
    print "z = ", dist.r[iterm, 2][flt]
    print "F  = ", MAXWELL_FACTOR * field[flt]
    print "E0 = ", external_field(dist.r[iterm, :][flt, :], t)
    print "SF = ", MAXWELL_FACTOR * sfields[flt]

    absE = sqrt(sum(E**2, axis=1))

    # An unit vector with the same direction as E
    u = E / absE[:, newaxis]

    # Now we can calculate the absolute value of the velocity
    if not SPRITES:
        vabs = TIP_MOBILITY * where(absE > TIP_MIN_FIELD, 
                                    absE - TIP_MIN_FIELD, 0)
    else:
        theta = sprite_theta(dist.r[iterm, :])
        vabs = TIP_MOBILITY * where(absE > TIP_MIN_FIELD * theta, 
                                    absE - TIP_MIN_FIELD * theta, 0) / theta
    v = u * vabs[:, newaxis]
    
    return v


def self_fields(tr, dist):
    """ Calculates the fields created by the charges at the streamer tips
    on themselves. """
    
    iterm = tr.terminals()
    parents = tr.parents()[iterm]
    
    dr = dist.r[iterm, :] - dist.r[parents, :]
    u = dr / (sqrt(sum(dr**2, axis=1)))[:, newaxis]
    
    theta2 = 1 if not SPRITES else sprite_theta(dist.r[iterm, :])[:, newaxis]**2

    # We take here the abs of q because the self-fields point always
    # to the advacement of the streamer.
    return (abs(dist.q[iterm][:, newaxis]) * u 
            / (CONDUCTOR_THICKNESS**2 / theta2))


def relax(box, tr, dist, t, dt):
    """ Relax the conductor :class:`tree.Tree` *tr* for a time *dt*. 

    Arguments:

      * *tr*: the :class:`tree.Tree` instance containing the tree structure.
      * *dist*: The tree distribution in a Distribution namedtuple.
      * *t*: Initial time of the time-step.
      * *dt*: the time step.
    """
    global latest_phi, error, error_dq
    
    #with ContextTimer("re-computing Ohm matrix"):
    # If we have an electrode, we fix q[0] by setting the first row of
    # M to zero.  
    fix = [] if ELECTRODE is None else [0]
    
    if SPRITES:
        factor = (sprite_theta(tr.midpoints(dist)) 
                  ** CONDUCTIVITY_SCALE_EXPONENT)
    else:
        factor = None

    # On Fri Aug 31 11:46:47 2012 I found a factor 2 here that I do not know
    # where it comes from.  Probably it was a reliq of the mid-points approach
    # (But it was duplicated in ohm_matrix anyway!).  I am removing it.
    M = CONDUCTANCE * tr.ohm_matrix(dist, factor=factor, fix=fix)
    n = len(dist.q)
    
    def f(t0, q):
        global latest_phi, error, error_dq

        phi = solve_phi(tree.Distribution(dist.r, q), box)
        # err = sqrt(sum((phi - box.phi)**2)) / len(phi)        
        latest_phi = phi
        error = phi - latest_phi
        error_dq = M.dot(error)
        
        dq = M.dot(phi + external_potential(dist.r, t + t0))

        return dq

    d = ode(f).set_integrator('vode',  nsteps=250000, rtol=1e-8)
    d.set_initial_value(dist.q, 0.0)
    d.integrate(dt)

    return tree.Distribution(r=dist.r, q=d.y)


def solve_phi(dist, box=None):
    if len(dist.q) >= FMM_THRESHOLD and box is not None:
        # with ContextTimer("FMM") as ct_fmm: 
        box.update_charges(dist.q)

        box.upward(MULTIPOLAR_TERMS)
        box.downward()
        box.solve_all(a=CONDUCTOR_THICKNESS, field=False)

        box.collect_solutions(field=False)

        phi = MAXWELL_FACTOR * box.phi
    else:
        # with ContextTimer("direct") as ct_direct:
        if ELECTRODE is None:
            distx = dist
        else:
            rx, qx = ELECTRODE.extend(dist.r, dist.q)
            distx = tree.Distribution(r=rx, q=qx)

        if not SPRITES:
            phi0 = MAXWELL_FACTOR * mpolar.direct(distx.r, distx.q, dist.r,
                                                  CONDUCTOR_THICKNESS)
        else:
            theta = sprite_theta(distx.r)
            phi0 = MAXWELL_FACTOR * mpolar.direct_ap(distx.r, distx.q, dist.r,
                                                CONDUCTOR_THICKNESS / theta)
        phi = phi0

    return phi


def symmetric_gaussian(dr, sigma):
    """ Samples a branch from a symmetric, gaussian branching model.
    In a plane perpendicular to dr we sample dr1 from a cylindrically
    symmetric gaussian distribution; the two branching points are dr1 and
    its symmetric vector wrt dr. """

    e1, e2 = perp_basis(dr)

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


def perp_basis(v):
    """ Returns two vectors e1, e2 such that (e1, e2, v/|v|) is an orthonormal
    basis. """

    u = v / norm(v)
    # We find two unit vectors orthonormal to u (also dr); note that this
    # fails if u is parallel to x !!!
    
    ex = array([1.0, 0, 0])

    e1 = ex - dot(u, ex) * u
    e1 = e1 / norm(e1)
    
    e2 = cross(u, e1)

    return e1, e2



def external_field(r, t):
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
    if SPRITES:
        return sprite_field(r, t)

    field = EXTERNAL_FIELD_VECTOR[newaxis, :]
    
    if ELECTRODE_POTENTIAL == 0:
        return field

    center = array([0.0, 0.0, ELECTRODE_RADIUS])
    dr = r - center[newaxis, :]
    rabs = sqrt(sum(dr**2, axis=1))
    field = field + (ELECTRODE_RADIUS * ELECTRODE_POTENTIAL
                     * dr / rabs[:, newaxis]**3)

    return field


def external_potential(r, t):
    """ Calculates the external potential at points *r*.  See above, in
    external_field for the risks here.
    """
    if SPRITES:
        return sprite_potential(r, t)

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

#
#  SIMULATION OF SPRITE DISCHARGES:
#  --------------------------------
#  The simulation of sprites differs from 'standard' simulations in two
#  aspects:
#  1. We allow an inhomogeneous density and rescale the relevant
#     streamer parameters accordingly
#  2. As external field, we use the field created by a charge sitting
#     between parallel electrodes.
#
#  Here we define the functions needed to implement this functionality.

def sprite_theta(r):
    """ Calculates the density ratio (theta=n(z)/n(0)) according
    to ``SCALE_HEIGHT``. """
    return exp(-r[:, 2] / SCALE_HEIGHT)


def sprite_instant_charge(t):
    """ Returns the cloud charge in a CG discharge at time t. """
    return (CHARGE_TOTAL / (TAU_1 - TAU_2) 
            * (+ TAU_2 * (exp(-t / TAU_2) - 1) 
               - TAU_1 * (exp(-t / TAU_1) - 1)))


def sprite_potential(r, t):
    """ Calculates the potential created at *r* by a charge *qt* inside a 
    parallel-plates capacitor. """
    qt = sprite_instant_charge(t)
    return MAXWELL_FACTOR * mpolar.direct(SPRITE_RS, qt * SPRITE_QS, r, 0.0)


def sprite_field(r, t):
    """ Calculates the field created at *r* by a charge *qt* inside a 
    parallel-plates capacitor. """
    qt = sprite_instant_charge(t)
    return MAXWELL_FACTOR * mpolar.field_direct(SPRITE_RS, 
                                                qt * SPRITE_QS, r, 0.0)


def sprite_charge_sources():
    """ Calculates the source charges of the applied field for sprite 
    simulations.  These are the mirror charges of the charge in the cloud
    created by reflections at the ground and the ionosphere. 
    We use mirror charges and a total of ``4*CHARGE_REFLECTIONS`` reflections.
    Returns *rs*, *qs*, where *rs* contains the charge locations
    and *qs* the charge sign (i.e. you have to multiply *qs* by
    q(t) to do anything meaningful. """

    # Location of the source charges
    # Source charge
    z0 = -(IONOSPHERE_HEIGHT - CHARGE_HEIGHT)
    
    # First reflection
    z1 = -z0

    # Each +/- charge repeats with a period 2 * IONOSPHERE_HEIGHT
    zpos = z0 + (2 * arange(-CHARGE_REFLECTIONS, CHARGE_REFLECTIONS + 1) 
                 * IONOSPHERE_HEIGHT)
    zneg = z1 + (2 * arange(-CHARGE_REFLECTIONS, CHARGE_REFLECTIONS + 1) 
                 * IONOSPHERE_HEIGHT)
    z = r_[zpos, zneg]

    qs = r_[ones((len(zpos),)), -ones((len(zneg),))]
    rs = zeros((len(z), 3))
    rs[:, 2] = z
    
    return rs, qs




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

