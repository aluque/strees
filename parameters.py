""" This is the parameter description file.  Basically it provides a namespace
where store a function per input parameter.  Each function receives a string
and is responsible for converting to the appropriate type, checking for
allowed values.  The docstring of the function is a description of the
parameter. """

from readinput import positive, nonnegative, default

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


@positive
def maxwell_factor(s):
    """ Maxwell factor for the potential and electric fields.  In SI units
    it is 1 / 4 pi epsilon_0"""
    return float(s)

    
@positive
def tip_mobility(s):
    """ Ratio between the tip velocity of each streamer and the local field. """
    return float(s)


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


@nonnegative
def branching_probability(s):
    """ Probability that a filament branches per unit distance travelled. """
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
def single_branching_z(s):
    """ If nonzero, performs a single branching at the given z. """
    return float(s)


@default(False)
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
