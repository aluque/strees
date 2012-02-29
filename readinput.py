""" Read an .ini input file.
    We use python's ConfigParser.
"""

import os, os.path, socket
from ConfigParser import SafeConfigParser
from warnings import warn

def guess_type(s):
    """ Converts s into int or float if possible.  If not, leave it as a
    string. This will give us problems if the user wants to e.g. use
    filenames or run names that can be parser as a number.  So have to
    implement a better approach whereby parameter names are associated with
    types.
    """

    try:
        return int(s)
    except ValueError:
        pass

    try:
        return float(s)
    except ValueError:
        pass

    return s


# These are decorators to check allowed values.
def positive(func):
    def f(s):
        r = func(s)
        if not r >= 0:
            raise ValueError("%s must be positive" % func.func_name)
        return r
    
    f.__doc__ = func.__doc__
    f.func_name = func.func_name

    return f


# These are decorators to check allowed values.
def nonnegative(func):
    def f(s):
        r = func(s)
        if r < 0:
            raise ValueError("%s must be positive" % func.func_name)
        return r
    
    f.__doc__ = func.__doc__
    f.func_name = func.func_name

    return f



def load_input(fname, parameters, d=None, upper=False):
    """ Loads an input file and stores its values in the dictionary
    d.   If upper is true, transforms the parameter names to upper case. """

    if d is None:
        d = dict()
    
    config = SafeConfigParser()
    defaults = dict(home=os.environ['HOME'],
                    user=os.getlogin(),
                    cwd=os.getcwd(),
                    input=os.path.splitext(os.path.basename(fname))[0]
                    hostname=socket.gethostname())
                    
    config.read(fname)

    for section in ['global', 'parameters']:
        for name, value in config.items(section, vars=defaults):
            try:
                func = getattr(parameters, name)
                value = func(value)
                print "%-30s =\t%-10s [%s]" % (name, repr(value), func.__doc__)
                
            except AttributeError:
                warn("'%s' is not defined as a parameter." % name)
                value = guess_type(value)

            if upper:
                name = name.upper()

            d[name] = value

    return d


def main():
    import sys

    d = load_input(sys.argv[1], d=globals(), upper=True)
    print TIME_STEP
    

if __name__ == '__main__':
    main()
    
