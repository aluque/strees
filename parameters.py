""" This is the parameter description file.  Basically it provides a namespace
where store a function per input parameter.  Each function receives a string
and is responsible for converting to the appropriate type, checking for
allowed values.  The docstring of the function is a description of the
parameter. """

from readinput import positive, nonnegative

def run_name(s):
    """ Name of the run. """
    return s

def out_file(s):
    """ File name (including path) of the .h5 output file. """
    return s

def desc(s):
    """ A comment to describe this simulation.  It is ignored by the code. """
    return s


@positive
def random_seed(s):
    """ Seed for the random generator. """
    return int(s)
    
    
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
    """ Probability per unit time that a filament branches. """
    return float(s)


@positive
def branching_sigma(s):
    """ Standard deviation of the branching dispacement in the symmetric gaussian branching model. """
    return float(s)


