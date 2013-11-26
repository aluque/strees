""" This is the parameter description file.  Basically it provides a namespace
where we store a function per input parameter.  Each function receives a string
and is responsible for converting to the appropriate type, checking for
allowed values.  The docstring of the function is a description of the
parameter. """

from readinput import positive, nonnegative, default, boolean

def run_name(s):
    """ Name of the run. """
    return s

def out_file(s):
    """ File name (including path) of the .h5 output file. """
    return s

def desc(s):
    """ A comment to describe this simulation.  It is ignored by the code. """
    return s


@default(-1)
def random_seed(s):
    """ Seed for the random generator. If < 0 the system time is used (default)."""
    return int(s)
    
def dummy(s):
    """ This parameter is completely ignored.  It is used only to produce
    series of identical runs. """

@positive
def max_charges_per_box(s):
    """ Maximum number of charges per box in the FMM refinement scheme. """
    return int(s)


@positive
def fmm_threshold(s):
    """ Threshold in the number of charges between using the direct solver \
    and FMM solver. """

    return int(s)


@positive
def multipolar_terms(s):
    """ Order of the multipole expansion in the FMM. """
    return int(s)


def external_field(s):
    """ Externally applied electric field in the z direction. """
    return float(s)


@positive
def conductance(s):
    """ Conductance of the channels. """
    return float(s)

@default(0)
@positive
def negative_conductance(s):
    """ Conductance of the negative channels. """
    return float(s)


@positive
def maxwell_factor(s):
    """ Maxwell factor for the potential and electric fields.  In SI units
    it is 1 / 4 pi epsilon_0"""
    return float(s)

    
@positive
def tip_mobility(s):
    """ Ratio between the tip velocity of each streamer and the local field. """
    return float(s)

@default(0)
@positive
def negative_tip_mobility(s):
    """ Ratio between the tip velocity of each streamer and the local field
    for negative streamers."""
    return float(s)

@default(0)
@nonnegative
def tip_min_field(s):
    """ Minimum field at the tip for a streamer to propagate. """
    return float(s)

@default(0)
@nonnegative
def initial_nodes(s):
    """ Starts the simulation with a vertical string with this number of 
    charged nodes separated a distance CONDUCTOR_THICKNESS. """
    return int(s)


@positive
def end_time(s):
    """ Final time of the simulation. """
    return float(s)

@positive
def time_step(s):
    """ Timestep of the simulation. """
    return float(s)


@positive
def conductor_thickness(s):
    """ Thickness of the conductors for the thin-wire approximation. """
    return float(s)


@default(0)
@positive
def negative_conductor_thickness(s):
    """ Thickness of the negative streamer channels. """
    return float(s)


@nonnegative
def branching_probability(s):
    """ Probability that a filament branches per unit distance travelled. """
    return float(s)


@default(1)
@nonnegative
def branching_radius_ratio(s):
    """ After branching, ratio of the descendant's radius relative to the parent """
    return float(s)

@default(1)
@nonnegative
def branching_conductance_ratio(s):
    """ After branching, ratio of the descendant's conductance relative to the parent """
    return float(s)


@positive
def branching_sigma(s):
    """ Standard deviation of the branching dispacement in the symmetric gaussian branching model. """
    return float(s)


@default(0)
@nonnegative
def single_branching_time(s):
    """ If nonzero, performs a single branching at the given time. """
    return float(s)


@default(0)
def threshold_z(s):
    """ If nonzero, changes the parameters reading from the [lower] section
    when the streamer reaches this altitude. """
    return float(s)

@default(0)
def fixed_branching_angle(s):
    """ If nonzero, fixes the angle between sibling branches. """
    return float(s)

@default(False)
@boolean
def branch_in_xz(s):
    """ If true, branches always within the XZ plane. """
    return (s.lower() == 'true')

@default('null')
def electrode_geometry(s):
    """ The electrode geometry. """
    return s.lower()


@positive
def electrode_radius(s):
    """ Radius of an spherical electrode. """
    return float(s)

@default(0)
def electrode_potential(s):
    """ Electrostatic potential of a (spherical) electrode. """
    return float(s)


@default(1e10)
def max_step(s):
    """ Longest step for a channel in dt.  The timestep will be reduced
    to satisfy this constraint"""
    return float(s)


@default(True)
@boolean
def end_with_reconnection(s):
    """ If true, finishes when a reconnection is detected """
    return (s.lower() == 'true')


#
# Parameters for the simulation of sprites
#
@default(False)
@boolean
def sprites(s):
    """ If true, implements the sprite-simulation functionality. """
    return (s.lower() == 'true')

@default(1e10)
def scale_height(s):
    """ Scale height of the densities. """
    return float(s)

@default(-1)
def conductivity_scale_exponent(s):
    """ Exponent of the scaling of streamer radius with n. """
    return float(s)


@default(1e-3)
@nonnegative
def trans_break_time(s):
    """ The typical time of transversal breakdown probability density. """
    return float(s)

@default(1e3)
@nonnegative
def trans_break_length(s):
    """ The typical length of transversal breakdown probability density. """
    return float(s)

@default(4.0)
@positive
def trans_break_alpha(s):
    """ Exponent of the transversal breakdown probability. """
    return float(s)

@default(0.0)
@nonnegative
def trans_break_field(s):
    """ Exponent of the transversal breakdown probability. """
    return float(s)


@default(0.0)
@nonnegative
def tau_1(s):
    """ Decay time scale of the cloud-to-ground current. """
    return float(s)

@default(0.0)
@nonnegative
def tau_2(s):
    """ Rise time scale of the cloud-to-ground current. """
    return float(s)

@default(0.0)
def ionosphere_height(s):
    """ Altitude of the ionosphere (upper electrode). """
    return float(s)

@default(0.0)
def charge_height(s):
    """ Altitude of the cloud charge. """
    return float(s)

@default(0.0)
def charge_total(s):
    """ Total charge transferred in the CG stroke (<0 for a +CG, >0 for a -CG). """
    return float(s)

@default(2)
@nonnegative
def charge_reflections(s):
    """ Number of charge reflection of each sign in each direction. """
    return int(s)

@default('')
def effective_ionization_rate_file(s):
    """ A file with the effective ionization rate, with columns E ~ nu. """
    return s
